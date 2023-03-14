import json, re, sys
from collections import defaultdict

def print_graphviz(f, d, l):
    # build dict of virtual "point" mapped to real hit ptr
    point2data = {}
    for item in l:  # collection of (virtual) "points" in iteration
        for vPoint in item:
            point2data[vPoint["ptr"]] = vPoint["c"]
    assert not (set(point2data.values()) ^ set(hitsDict.keys())) # -> empty set
    # build dict of links
    links = dict()
    for item in l:  # collection of (virtual) "points" in iteration
        for vPoint in item:
            for link in vPoint['refs']:
                kFrom = vPoint['c']
                kTo = point2data[link['link'][0]]
                state = link['state']
                links[(kFrom, kTo)] = state
    # ---
    f.write('digraph G {\n')  # "    splines = false;\n"
    byLayers = defaultdict(list)
    # put basic connections
    for ptr, item in d.items():
        rank = int(item['layer'])
        byLayers[rank].append(ptr)
        #for con in item['cons']:
        #    assert (ptr, con) in links.keys()
        #    #print('link: ', (ptr, con) )
        #    f.write('    \"%s\" -> \"%s\" [label=\"%d\"];\n'%(ptr, con, links[(ptr, con)]))
        for con, st in ([linkPair, linkState] for linkPair, linkState in links.items() if linkPair[0] == ptr):
            if st == 1: continue
            f.write('    \"%s\" -> \"%s\" [label=\"%d\"];\n'%(ptr, con[1], st))
    # specify rank
    for layerNo, items in byLayers.items():
        f.write('    { rank = same ; %s}\n'%(', '.join(f'"{item}"' for item in items)) )
    # node attributes
    for ptr, item in d.items():
        f.write('   \"%s\" [label=\"%s\"]\n'%(ptr, item['label']))
    f.write('}\n')

#with open('dict.json') as f:
#    hitsDict = json.load(f)

with open(sys.argv[1]) as f:
    dbgLog = json.load(f)

# split labels to get layer and hit ID
rxsLbl = re.compile('^\\((?P<layer>\d+),(?P<hitNo>\d+)\\)')
hitsDict = {}
for ptr, item in dbgLog['dict'].items():
    m = rxsLbl.match(item)
    #item.update(m.groupdict())
    hitsDict[ptr] = dict(m.groupdict())
    hitsDict[ptr]['label'] = item

print_graphviz(sys.stdout, hitsDict, dbgLog['its'][-1])
#print(dbgLog['iterations'][0])  # XXX
