import sys
import re

verbal = True
version = ''
warns = {}

def varsub(ys, v=True):
    global verbal
    global version
    verbal = v
    version = sys.version[0]

    hit = 1

    while hit:
        hit = 0
        objs = [ys]

        while objs:
            obj = objs.pop()
            t = str(type(obj))
            if isinstance(obj, dict):
                keys = obj.keys()

                for key in keys:
                    t = str(type(obj[key]))
                    if isinstance(obj[key], str):
                        if __scan4vars(ys, obj, key):
                            hit = 1
                    elif isinstance(obj[key], list):
                        objs.append(obj[key])
                    elif isinstance(obj[key], dict):
                        objs.append(obj[key])
            else:  # must be list
                for idx in range(len(obj)):
                    el = obj[idx]
                    if isinstance(el, str):
                        if __scan4vars(ys, obj, idx):
                            hit = 1
                    elif isinstance(el, list):
                        objs.append(el)
                    elif isinstance(el, dict):
                        objs.append(el)
    if verbal and warns:
        for k in warns:
            print(warns[k], file=sys.stderr)

def __scan4vars(ys, obj, key):
    pat = re.compile(r'\${?([\w\d]+)}?')
    full = re.compile(r'\${?([\w\d]+)}?$')
    ps = {}
    hit = 0
    for m in re.findall(pat, obj[key]):
        ps[m] = 1  # handle duplicated matches
    if ps:
        if re.match(full, obj[key]):  # not embedded substitution
            k = list(ps.keys())[0]
            if k in ys:
                obj[key] = ys[k]
                hit = 1
            elif verbal:  # warning, quitting?
                warns[k] = k + ' in "%s" not defined' % (obj[key])
        else:
            for k in ps.keys():
                if k in ys:
                    p = r'\${?%s}?' % k
                    if isinstance(ys[k], str):
                        obj[key] = re.sub(p, ys[k], obj[key])
                        hit = 1
                    elif isinstance(ys[k], (int, float)):
                        obj[key] = re.sub(p, str(ys[k]), obj[key])
                        hit = 1
                    elif version == 2 and isinstance(ys[k], long):
                        obj[key] = re.sub(p, str(ys[k]), obj[key])
                        hit = 1
                    elif verbal:
                        warns[k] = k + ' is not a valid type for substitution in "%s"' % (obj[key])
                elif verbal:  # warning, quitting?
                    warns[k] = k + ' in "%s" not defined' % (obj[key])
    return hit
