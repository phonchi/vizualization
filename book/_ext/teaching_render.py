"""Keep identical inline assets from rendering twice with Jupyter Book 1.x.

Sphinx can register a late inline asset with filename=None, then replay the
registry with filename="".  Those assets render identically but compare unequal.
Use the public page-context hook to normalize only this duplicate case; leave
external scripts and inline scripts with different attributes untouched.
"""

from typing import Any


def deduplicate_inline_scripts(
    app: Any,
    pagename: str,
    templatename: str,
    context: dict[str, Any],
    doctree: Any,
) -> None:
    """Remove exact inline duplicates without changing their first-use order."""
    scripts = context.get("script_files")
    if scripts is None:
        return

    seen = set()
    unique = []
    for script in scripts:
        filename = getattr(script, "filename", None)
        attributes = getattr(script, "attributes", None)
        # Unknown template objects and linked assets are outside this fix.
        if filename not in (None, "") or not isinstance(attributes, dict):
            unique.append(script)
            continue
        if "body" not in attributes:
            unique.append(script)
            continue

        # Include every attribute and priority, not only body: an inline module
        # and a classic script containing the same text are distinct assets.
        key = (
            filename or "",
            getattr(script, "priority", None),
            tuple(sorted((name, repr(value)) for name, value in attributes.items())),
        )
        if key not in seen:
            seen.add(key)
            unique.append(script)

    context["script_files"] = unique


def setup(app: Any) -> dict[str, Any]:
    # Run after normal theme context callbacks, while preserving Sphinx's final
    # priority sort and normal script-tag generation.
    app.connect("html-page-context", deduplicate_inline_scripts, priority=1000)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
