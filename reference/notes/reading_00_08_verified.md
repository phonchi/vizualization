# 首頁與第一部延伸閱讀核實記錄

核對日期：2026-09-09。章節依 `book/_toc.yml`，而非舊版 README 的章數。
以下區塊供教材章末使用；文獻導讀是依章節需求撰寫的閱讀建議。

## 00_intro

```markdown
## 參考資料與延伸閱讀

- **入門・免費網站**：USGS，〈[Can you predict earthquakes?](https://www.usgs.gov/faqs/can-you-predict-earthquakes)〉。
  先讀這篇短文，釐清「預測一次地震」與「估計地震機率」的差別，也想想為什麼事後找到異常還不能算預測成功。
- **台灣觀測背景・免費文章**：蕭乃祺（2019），〈[臺灣地震測報的發展](https://www.ntsec.edu.tw/liveSupply/detail.aspx?a=6829&cat=6841&lid=16154&p=1&print=1)〉，《科學研習》58 卷 6 期。
  從地震儀、GNSS 與地下水觀測的建置歷程，認識本站資料從何而來。文中的測站數與作業現況是當時的紀錄，可與本站測站清單對照。
- **進階・開放取用論文**：Mizrahi et al.（2024），〈[Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823)〉，*Reviews of Geophysics*。
  建議先讀 Plain Language Summary，再看模型開發、檢驗與溝通三部分的關係；這正是第二部各章串起來要回答的問題。
```

## 01_overview

```markdown
## 參考資料與延伸閱讀

- **入門・免費文章**：蕭乃祺（2019），〈[臺灣地震測報的發展](https://www.ntsec.edu.tw/liveSupply/detail.aspx?a=6829&cat=6841&lid=16154&p=1&print=1)〉，《科學研習》58 卷 6 期。
  優先讀「基礎」中的觀測網介紹，對照本章地圖上的地震、GNSS 與地下水測站；留意歷史測站數不等於目前仍在運作的站數。
- **資料來源・免費網站**：中央氣象署／FDSN，〈[Central Weather Administration Seismographic Network（T5）](https://www.fdsn.org/networks/detail/T5/)〉。
  查閱地震測網的代碼、營運機構與資料引用方式，理解測網代碼與單一測站代碼的差別。這是地震測網的登錄資料，不代表本站全部地球物理測網。
- **系統背景・免費官方文章**：中央氣象署地震測報中心，〈[資訊服務－地球物理資料管理系統](https://scweb.cwa.gov.tw/zh-tw/page/twenty/129)〉，《20 周年專刊》。
  說明不同觀測系統的資料如何整合、保存與供人查詢，可接著思考本章的測網沿革為何會影響資料涵蓋範圍。
```

## 02_download

```markdown
## 參考資料與延伸閱讀

- **操作入口・免費網站**：中央氣象署，〈[臺灣地震與地球物理資料管理系統（GDMS）](https://gdms.cwa.gov.tw/)〉。
  對照本章的查詢、申請與取件流程，查看官方資料項目及使用說明；資料下載需先註冊登入。
- **格式查詢・免費文件**：FDSN，〈[Channel codes](https://docs.fdsn.org/projects/source-identifiers/en/v1.0/channel-codes.html)〉。
  查閱 band、source 與 subsource codes，對照本章的 HHZ、HH1 與 HH2。重點是頻道代碼描述什麼訊號，不要只把它當成檔名的一部分。
- **核心論文**：Beyreuther et al.（2010），〈[ObsPy: A Python Toolbox for Seismology](https://doi.org/10.1785/gssrl.81.3.530)〉，*Seismological Research Letters*；[出版學會免費全文](https://www.seismosoc.org/Publications/SRL/SRL_81/srl_81-3_es/)。
  了解下載後的波形為什麼適合交給 ObsPy 統一讀取與處理。論文中的程式介面屬於早期版本，實際寫程式時搭配[目前的官方教學](https://docs.obspy.org/tutorial/)。
```

## 03_groundwater

```markdown
## 參考資料與延伸閱讀

- **入門・免費網站**：USGS，〈[How does an earthquake affect groundwater levels and water quality in wells?](https://www.usgs.gov/faqs/how-does-earthquake-affect-groundwater-levels-and-water-quality-wells)〉。
  先分清震波經過時的水位振盪、震後水位偏移與水質變化，再回頭判斷本章各口井究竟觀測到了哪一種反應。
- **觀測機制・免費官方教材**：USGS（2003），〈[Earthquakes—Rattling the Earth's Plumbing System](https://pubs.usgs.gov/fs/fs-096-03/)〉，Fact Sheet 096-03。
  用不同地點的觀測說明地震如何擾動地下水，適合延伸本章「同一場地震，不同井反應不同」的討論。
- **核心綜述・全文需訂閱**：Roeloffs, E. A.（1988），〈[Hydrologic precursors to earthquakes: A review](https://doi.org/10.1007/BF00878996)〉，*Pure and Applied Geophysics*。
  重點讀潮汐響應、承壓含水層與應變的關係，以及如何排除氣壓、降雨與抽水影響；文獻整理的候選前兆不能直接當成已驗證的預測方法。
```

## 04_geomagnetic

```markdown
## 參考資料與延伸閱讀

- **背景查詢・免費網站**：GFZ，〈[Kp index](https://kp.gfz.de/en/)〉。
  查閱全球地磁活動指數與說明，對照本章的擾動日；先檢查太空天氣背景，再討論單站或兩站差值中的局部異常。
- **入門・免費官方文章**：USGS Geomagnetism Program，〈[Overview](https://www.usgs.gov/programs/geomagnetism/science/overview)〉。
  認識地磁觀測的研究用途，並閱讀其中對地震相關磁場宣稱的檢驗，延伸本章對參考站、背景場與儀器干擾的討論。
- **核心綜述・全文可能需訂閱**：Johnston, M. J. S.（1997），〈[Review of electric and magnetic fields accompanying seismic and volcanic activity](https://doi.org/10.1023/A:1006500408086)〉，*Surveys in Geophysics*；[USGS 免費摘要](https://www.usgs.gov/publications/review-electric-and-magnetic-fields-accompanying-seismic-and-volcanic-activity)。
  對照壓磁效應、流體相關電磁效應與觀測頻帶，特別區分同震訊號和震前訊號；看到地震伴隨的磁場變化，不等於能據此提前預測地震。
```

## 05_seismic

```markdown
## 參考資料與延伸閱讀

- **操作入門・免費教材**：ObsPy 開發團隊，〈[Tutorial](https://docs.obspy.org/tutorial/)〉。
  依序看 Reading Seismograms、Filtering Seismograms 與 Plotting Spectrograms，把本章的讀檔、濾波及時頻圖流程接起來。
- **工具論文**：Beyreuther et al.（2010），〈[ObsPy: A Python Toolbox for Seismology](https://doi.org/10.1785/gssrl.81.3.530)〉，*Seismological Research Letters*；[出版學會免費全文](https://www.seismosoc.org/Publications/SRL/SRL_81/srl_81-3_es/)。
  了解 ObsPy 如何把不同格式的地震資料接到共同的處理流程；論文用來認識設計背景，程式語法以目前官方文件為準。
- **台灣規模統計・論文**：Wang et al.（2015），〈[b-Values Observations in Taiwan: A Review](https://doi.org/10.3319/TAO.2015.04.28.01%28T%29)〉，*Terrestrial, Atmospheric and Oceanic Sciences*。
  延伸本章的規模－次數圖與震前震後比較，留意目錄品質、取樣區間和構造差異如何影響 $b$ 值；不同的估計值本身還不足以證明前兆。
- **台灣餘震統計・論文**：Wang et al.（2016），〈[Studies on Aftershocks in Taiwan: A Review](https://doi.org/10.3319/TAO.2016.09.12.01)〉，*Terrestrial, Atmospheric and Oceanic Sciences*；[期刊文章與全文入口](https://tao.cgu.org.tw/index.php/articles/archive/geophysics/item/1498-2016091201t)。
  對照台灣不同序列的餘震分布、Omori 衰減與觸發機制，思考本章花蓮個案的結果能推廣到什麼範圍。
```

## 06_gnss

```markdown
## 參考資料與延伸閱讀

- **資料實作・免費網站**：Nevada Geodetic Laboratory，〈[Plug and Play GPS Data Products](https://geodesy.unr.edu/PlugNPlayPortal.php)〉。
  從測站清單進入位置時間序列與資料格式說明，延伸本章「原始衛星觀測」和「已解算座標」的區別；比較位移前先確認參考框架及單位。
- **格式查詢・免費文件**：IGS／RTCM RINEX Working Group，〈[RINEX](https://igs.org/wg/rinex/)〉。
  依手上檔案的版本選擇規格，查閱標頭與觀測量代碼，對照本章的 RINEX 解析；檔頭的近似座標不等於逐時刻解算的位置序列。
- **資料方法導讀・免費文章**：Blewitt, G., Hammond, W. C., & Kreemer, C.（2018），〈[Harnessing the GPS Data Explosion for Interdisciplinary Science](https://doi.org/10.1029/2018EO104623)〉，*Eos*；[免費全文](https://eos.org/science-updates/harnessing-the-gps-data-explosion-for-interdisciplinary-science)。
  了解大量 GNSS 資料如何整理成可研究的速度與位移產品，特別看季節變化、地震階變及設備更動如何影響時間序列的解讀。
```

## 07_case_hualien2024

```markdown
## 參考資料與延伸閱讀

- **案例背景・免費報告**：臺灣地震科學中心（TEC）（2024），〈[2024 M7.2 花蓮地震](https://tec.earth.sinica.edu.tw/specialEQ/pdf/20240403hwalien.pdf)〉。
  先看主震位置、區域構造與觀測摘要，再與本章的事件時間軸對照；這是事件報告，閱讀時要留意版本與資料來源。
- **對照資料・免費網站**：USGS（2024），〈[Hualien earthquake：事件頁 us7000m9g4](https://earthquake.usgs.gov/earthquakes/eventpage/us7000m9g4/executive)〉。
  比較官方事件參數、震源機制與地動產品；USGS 與本章使用的規模尺度、定位結果及時間表示可能不同，先核對定義再比較數字。
- **核心案例・開放取用論文**：Zheng et al.（2024），〈[Thrust-dominated unilateral rupture of a blind listric fault associated with the 2024 Hualien earthquake](https://doi.org/10.1038/s41598-024-82971-x)〉，*Scientific Reports*。
  看作者如何結合地震波與大地測量資料推估破裂過程，特別注意不同觀測對斷層幾何的約束；這是超越本章時間序列對照、進一步建立震源模型的例子。
- **台灣比較案例・論文**：Wang et al.（2016），〈[Studies on Aftershocks in Taiwan: A Review](https://doi.org/10.3319/TAO.2016.09.12.01)〉，*Terrestrial, Atmospheric and Oceanic Sciences*；[期刊文章與全文入口](https://tao.cgu.org.tw/index.php/articles/archive/geophysics/item/1498-2016091201t)。
  把花蓮序列放進台灣其他地震的背景中，對照餘震分布與觸發機制；這篇早於 2024 年，提供的是比較框架，並非花蓮事件的分析結果。
```

## 08_explore_ideas

```markdown
## 參考資料與延伸閱讀

- **入門・免費網站**：USGS，〈[Can you predict earthquakes?](https://www.usgs.gov/faqs/can-you-predict-earthquakes)〉。
  用短文檢查本章的核心問題：判準是否事先定好、沒有地震時是否也會出現相同異常，以及怎樣才算可檢驗的預測。
- **水文綜述・全文需訂閱**：Roeloffs, E. A.（1988），〈[Hydrologic precursors to earthquakes: A review](https://doi.org/10.1007/BF00878996)〉，*Pure and Applied Geophysics*。
  對照本章「背景扣乾淨了嗎」這一關，閱讀地下水異常的物理解釋與非構造因素，思考需要哪些額外觀測才能排除替代解釋。
- **電磁綜述・全文可能需訂閱**：Johnston, M. J. S.（1997），〈[Review of electric and magnetic fields accompanying seismic and volcanic activity](https://doi.org/10.1023/A:1006500408086)〉，*Surveys in Geophysics*；[USGS 免費摘要](https://www.usgs.gov/publications/review-electric-and-magnetic-fields-accompanying-seismic-and-volcanic-activity)。
  比較候選物理機制與實際觀測限制，延伸本章對儀器干擾、多站比較，以及同震效應不能當作震前證據的討論。
- **跨觀測綜述・全文可能需訂閱**：Cicerone, R. D., Ebel, J. E., & Britton, J.（2009），〈[A systematic compilation of earthquake precursors](https://doi.org/10.1016/j.tecto.2009.06.008)〉，*Tectonophysics*。
  看不同候選前兆如何整理成時間、空間、振幅與訊雜比等可比較欄位；彙整曾被報告的現象，和證明能前瞻預報，是兩個不同層次的工作。
- **統計方法・免費論文**：Wasserstein, R. L., & Lazar, N. A.（2016），〈[The ASA Statement on p-Values: Context, Process, and Purpose](https://doi.org/10.1080/00031305.2016.1154108)〉，*The American Statistician*。
  延伸本章對事後調整門檻與選擇性報告的提醒；重點看為什麼單一 $p$ 值不能取代研究設計、效果大小與完整的分析紀錄。
```

## 核實與原始素材

- 已盤點 `reference/`（含 Git 忽略的 PDF）、讀取 `notes/INDEX.md`、`taiwan.md`、`oef.md` 與各章內容。
- `Taiwan/臺灣預報發展.pdf`：全文核對蕭乃祺署名、三個測網、地下水/GNSS 建置沿革；與科教館原文一致。2019 年另由科教館第 58 卷 6 期目錄及氣象署出版品書目交叉核實。
- `Taiwan/[2015] b-Values Observations in Taiwan A Review.pdf`：首頁核對 Wang、Chen、Leu、Chang 四位作者、2015 年及 DOI；摘要和導論對應第 5 章。
- `Taiwan/[2016] Studies on Aftershocks in Taiwan A Review.pdf`：首頁核對作者、2016 年、DOI、摘要；第 5、7 章用其台灣餘震比較框架。期刊官方文章頁核實標題與摘要。
- `[2024] OEF_Review.pdf`：首頁及 Plain Language Summary 核對 Mizrahi et al. 2024，DOI 10.1029/2023RG000823；首頁導讀指向原文的三支柱。
- 第 8 章原有 Roeloffs 1988、Johnston 1997、Cicerone et al. 2009，已保留並補全正式書目及導讀。Roeloffs 由 Springer 原始頁核實；Johnston 由 USGS 作者機構紀錄核實；Cicerone 由 ScienceDirect 原始搜尋結果核實。
- 新補的資料操作來源由 GDMS、FDSN、ObsPy、IGS、NGL 官方頁核實；Beyreuther et al. 2010 另有出版學會原文，勿照抄早期 API。
- Zheng et al. 由 Scientific Reports 原始全文核實作者、2024-12-28 發表日、開放取用與方法；TEC 原始 PDF 提供中文事件導讀。
- 部分網站對自動擷取逾時或拒絕；「書目已核實」與「此環境的 HTTP 請求成功」分開記錄於驗證結果，不把搜尋結果當成全文閱讀。
- 以上不使用審稿中材料，也不複製或公開本機 PDF。
