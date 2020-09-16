"""A very simple function for fitting output to a console"""

def fit_to_console(obj, initial_indent='', subsequent_indent='', width=None,):
    """Return a string formatted to fit nicely in the console

    This is essentially a wrapper around textwrap and pprint; it pretty-prints the input object,
    then wraps it to the width of the output console.  The returned object is a string that can be
    printed without further adjustment.

    Parameters
    ----------
    obj: object
        Python object to be printed
    initial_indent: str [defaults to '']
        String that will be prepended to the first line of wrapped output.  Counts towards the
        line's width.  Note that this need not be whitespace; it could be a label, for example, to
        explain what `obj` is.
    subsequent_indent: str [defaults to '']
        String that will be prepended to all lines save the first of wrapped output; also counts
        towards each line's width.
    width: int or None [defaults to None]
        Full width of text to be output.  If None, the terminal size is detected, and the number of
        columns (minus 1) is used.

    """
    import shutil
    import textwrap
    import pprint
    columns, lines = shutil.get_terminal_size()
    full_width = width or columns-1
    wrapper = textwrap.TextWrapper(width=full_width, initial_indent=initial_indent, subsequent_indent=subsequent_indent)
    text = pprint.pformat(obj)
    return wrapper.fill(text)
