"""GDMS 會員登入 session。

GDMS 下載資料需要會員身分。本模組提供兩種方式建立已登入的 requests.Session：

1. 帳號密碼登入（.env 設 GDMS_USER / GDMS_PASS）——若登入頁啟用 reCAPTCHA
   可能失敗，此時改用方式 2。
2. 瀏覽器 cookie（.env 設 GDMS_COOKIE=PHPSESSID=xxxx）——先用瀏覽器登入
   GDMS，再從開發人員工具複製 cookie（步驟見教學第 2 章）。
"""

import os
import ssl
from pathlib import Path

import requests
from dotenv import load_dotenv
from requests.adapters import HTTPAdapter

BASE = "https://gdms.cwa.gov.tw"


class _GovCertAdapter(HTTPAdapter):
    """GDMS 的政府憑證缺 Subject Key Identifier 欄位，Python 3.13 預設的
    嚴格檢查（VERIFY_X509_STRICT）會拒絕連線。這裡仍驗證憑證鏈與主機名，
    只關閉嚴格模式。"""

    def init_poolmanager(self, *args, **kwargs):
        ctx = ssl.create_default_context()
        ctx.verify_flags &= ~ssl.VERIFY_X509_STRICT
        kwargs["ssl_context"] = ctx
        return super().init_poolmanager(*args, **kwargs)


def make_session() -> requests.Session:
    """建立可連 GDMS 的 requests.Session（處理憑證相容性）。"""
    s = requests.Session()
    s.headers["User-Agent"] = "gdms-toolkit-edu/1.0"
    s.mount("https://gdms.cwa.gov.tw", _GovCertAdapter())
    return s

# 從專案根目錄載入 .env
load_dotenv(Path(__file__).resolve().parent.parent / ".env")


class GDMSSession:
    """已登入的 GDMS session。用法：

    >>> gdms = GDMSSession()      # 自動讀 .env
    >>> gdms.logged_in
    True
    """

    def __init__(self, user: str | None = None, password: str | None = None,
                 cookie: str | None = None):
        self.session = make_session()
        self.logged_in = False

        cookie = cookie or os.getenv("GDMS_COOKIE")
        user = user or os.getenv("GDMS_USER")
        password = password or os.getenv("GDMS_PASS")

        if cookie:
            self._login_with_cookie(cookie)
        elif user and password:
            self._login_with_password(user, password)
        else:
            raise RuntimeError(
                "找不到 GDMS 帳號資訊。請在專案根目錄建立 .env，"
                "設定 GDMS_USER / GDMS_PASS（或 GDMS_COOKIE），"
                "範本見 .env.example。"
            )

    def _login_with_cookie(self, cookie: str):
        for pair in cookie.split(";"):
            name, _, value = pair.strip().partition("=")
            self.session.cookies.set(name, value, domain="gdms.cwa.gov.tw")
        self.logged_in = self._check()
        if not self.logged_in:
            raise RuntimeError("GDMS_COOKIE 無效或已過期，請重新從瀏覽器複製。")

    def _login_with_password(self, user: str, password: str):
        self.session.get(f"{BASE}/login.php", timeout=30)  # 先取得 session cookie
        r = self.session.post(
            f"{BASE}/php/loginProcess.php",
            data={"username": user, "password": password,
                  "g-recaptcha-response": "", "img-captcha": ""},
            timeout=30,
        )
        r.raise_for_status()
        try:
            self.logged_in = r.json().get("status") == 1
        except ValueError:
            self.logged_in = False
        if not self.logged_in:
            raise RuntimeError(
                f"帳密登入失敗（GDMS 回應：{r.text[:120]}）。"
                "請確認 .env 的 GDMS_USER / GDMS_PASS；若被驗證碼擋下，"
                "改用瀏覽器登入後複製 cookie 設定 GDMS_COOKIE，見教學第 2 章。"
            )

    def _check(self) -> bool:
        """以首頁是否出現登出鈕判斷 session 是否有效。"""
        r = self.session.get(f"{BASE}/index.php", timeout=30)
        return "logout.php" in r.text

    def get(self, url, **kw):
        return self.session.get(url, timeout=kw.pop("timeout", 60), **kw)

    def post(self, url, **kw):
        return self.session.post(url, timeout=kw.pop("timeout", 60), **kw)
