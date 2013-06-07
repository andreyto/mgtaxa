
def parseOrfOnContigIdApis(id):
    ret = {}
    ret["id_cont"],ret["start_orf"],ret["end_orf"],ignore = id.strip().split('_')
    ret["start_orf"] = int(ret["start_orf"])
    ret["end_orf"] = int(ret["end_orf"])
    if ret["end_orf"] < ret["start_orf"]:
        tmp = ret["end_orf"]
        ret["end_orf"] = ret["start_orf"]
        ret["start_orf"] = tmp
        ret["strand_orf"] = -1
    else:
        ret["strand_orf"] = 1
    ret["start_orf"] -= 1 # unit-offset is input, zero-offset is output
    return ret

