import re
import difflib
from mkdocs.plugins import BasePlugin


class AdmonitionTranslationPlugin(BasePlugin):
    """
    Convert GitHub-flavored admonitions to MkDocs numpy-style admonitions.
    """

    def on_page_markdown(self, markdown, page, config, files):
        lines = markdown.splitlines(keepends=True)
        out = []
        i = 0
        while i < len(lines):
            m = re.match(r'> \[\!(?P<tag>\w+)\]\s*$', lines[i])
            if m:
                tag = m.group('tag').lower()
                # emit the MkDocs admonition start
                out.append(f"!!! {tag}\n")
                i += 1
                # consume all the > … lines as the body
                while i < len(lines) and lines[i].startswith('>'):
                    # strip leading "> " (or just ">") and rstrip newline
                    body = lines[i].lstrip('> ').rstrip('\n')
                    out.append(f"    {body}\n")
                    i += 1
                continue
            # otherwise, copy line verbatim
            out.append(lines[i])
            i += 1

        return "".join(out)
