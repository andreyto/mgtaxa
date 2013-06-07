"""Methods for helping process JSON"""

def json_records_iterator(f,parse=True,parser=None,withJson=False):
    """Iterate over JSON records in file-like stream.
    @pre input is a sequence of valid JSON records
    separated by new line. Each record itself can be
    multiline. Each record is either array or object.
    @param f line iterator (e.g. file object)
    @return this generator method will yield JSON strings
    """
    if parse:
        if not parser:
            import json
            if withJson:
                parser = lambda s: (json.loads(s),s)
            else:
                parser = json.loads
    else:
        parser = lambda s: s
    buff = ""
    for s in f:
        s = s.strip()
        if buff and s and (buff[-1] in "]}" and s[0] in "[{"):
            yield parser(buff)
            buff = s
        else:
            buff += s
    if buff:
        yield parser(buff)

