"""Build a compact routable graph from the convoy network GeoJSON.

Output (little-endian binary, gzip):
  header: magic 'CVG1', nNodes u32, nEdges u32, nGeom u32
  nodes:  nNodes * (lon i32, lat i32)          coords * 1e5
  edges:  nEdges * (a u32, b u32, len f32, level u8, pad u8[3], gStart u32, gCount u32)
  geom:   nGeom  * (lon i32, lat i32)          interior points of each edge
Levels: 0=1TE, 1=2TE48, 2=TE72, 3=TE94, 4=TE120
"""
import json, math, struct, sys, gzip
from collections import defaultdict

SRC, OUT = sys.argv[1], sys.argv[2]
LEVEL = {'1TE': 0, '2TE48': 1, 'TE72': 2, 'TE94': 3, 'TE120': 4}
SNAP = 1e4        # node identity grid (~10 m)
BRIDGE_M = 25     # connect dead ends to a nearby node within this distance
SIMPLIFY_M = 8    # Douglas-Peucker tolerance

d = json.load(open(SRC))
key_of = {}   # snapped key -> node id
coords = []   # node id -> (lon, lat) float

def node(lon, lat):
    k = (round(lon * SNAP), round(lat * SNAP))
    i = key_of.get(k)
    if i is None:
        i = len(coords); key_of[k] = i; coords.append((lon, lat))
    return i

seg = {}  # (a,b) a<b -> level (max)
for f in d['features']:
    lvl = LEVEL.get(f['properties'].get('reseau_te'))
    if lvl is None:
        continue
    for line in f['geometry']['coordinates']:
        prev = None
        for lon, lat in line:
            n = node(lon, lat)
            if prev is not None and prev != n:
                k = (prev, n) if prev < n else (n, prev)
                if seg.get(k, -1) < lvl:
                    seg[k] = lvl
            prev = n
print('raw nodes', len(coords), 'raw segments', len(seg))

def dist(i, j):
    (x1, y1), (x2, y2) = coords[i], coords[j]
    kx = math.cos(math.radians((y1 + y2) / 2)) * 111320
    return math.hypot((x2 - x1) * kx, (y2 - y1) * 110540)

adj = defaultdict(dict)
for (a, b), l in seg.items():
    adj[a][b] = l; adj[b][a] = l

# --- bridge dead ends to nearby nodes (small digitising gaps)
grid = defaultdict(list)
G = 0.0005
for i in adj:
    x, y = coords[i]
    grid[(int(x // G), int(y // G))].append(i)
bridged = 0
for i in [n for n in adj if len(adj[n]) == 1]:
    x, y = coords[i]
    gx, gy = int(x // G), int(y // G)
    best, bd = None, BRIDGE_M
    nb = next(iter(adj[i]))
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in grid[(gx + dx, gy + dy)]:
                if j == i or j == nb or j in adj[i]:
                    continue
                dd = dist(i, j)
                if dd < bd:
                    best, bd = j, dd
    if best is not None:
        l = min(adj[i][nb], max(adj[best].values()))
        adj[i][best] = l; adj[best][i] = l; bridged += 1
print('bridged dead ends', bridged)

# --- contract degree-2 chains with uniform level
def is_junction(n):
    nb = adj[n]
    return len(nb) != 2 or len(set(nb.values())) != 1

def simplify(pts):
    if len(pts) < 3:
        return pts
    keep = [False] * len(pts); keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        s, e = stack.pop()
        (x1, y1), (x2, y2) = pts[s], pts[e]
        kx = math.cos(math.radians(y1)) * 111320; ky = 110540
        ax, ay = x1 * kx, y1 * ky; bx, by = x2 * kx, y2 * ky
        vx, vy = bx - ax, by - ay; L2 = vx * vx + vy * vy
        best, bi = 0, -1
        for k in range(s + 1, e):
            px, py = pts[k][0] * kx, pts[k][1] * ky
            if L2 == 0:
                dd = math.hypot(px - ax, py - ay)
            else:
                t = max(0, min(1, ((px - ax) * vx + (py - ay) * vy) / L2))
                dd = math.hypot(px - ax - t * vx, py - ay - t * vy)
            if dd > best:
                best, bi = dd, k
        if best > SIMPLIFY_M:
            keep[bi] = True; stack += [(s, bi), (bi, e)]
    return [p for p, k in zip(pts, keep) if k]

junction_ids = {}
def jid(n):
    if n not in junction_ids:
        junction_ids[n] = len(junction_ids)
    return junction_ids[n]

visited = set()  # undirected raw segments
edges = []       # (a, b, length, level, geom list)

def walk(start, nxt):
    lvl = adj[start][nxt]
    path = [start, nxt]
    visited.add((min(start, nxt), max(start, nxt)))
    prev, cur = start, nxt
    while not is_junction(cur) and cur != start:
        nb = [x for x in adj[cur] if x != prev]
        if not nb:
            break
        prev, cur = cur, nb[0]
        sk = (min(prev, cur), max(prev, cur))
        if sk in visited:
            break
        visited.add(sk)
        path.append(cur)
    length = sum(dist(path[k], path[k + 1]) for k in range(len(path) - 1))
    return path, length, lvl

for n in list(adj):
    if not is_junction(n):
        continue
    for m in list(adj[n]):
        if (min(n, m), max(n, m)) in visited:
            continue
        path, length, lvl = walk(n, m)
        geom = simplify([coords[p] for p in path])[1:-1]
        edges.append((jid(path[0]), jid(path[-1]), length, lvl, geom))
# isolated loops with no junction
for n in list(adj):
    for m in list(adj[n]):
        if (min(n, m), max(n, m)) not in visited:
            path, length, lvl = walk(n, m)
            geom = simplify([coords[p] for p in path])[1:-1]
            edges.append((jid(path[0]), jid(path[-1]), length, lvl, geom))

nodes = [None] * len(junction_ids)
for raw, j in junction_ids.items():
    nodes[j] = coords[raw]
edges = [e for e in edges if e[0] != e[1] or e[2] > 0]
# --- gap connectors: per vehicle class, join a dead end to the nearest node of
# another connected piece within GAP_M. Flag bit 3 marks them (off-network link).
GAP_M = 1500
GAP_FLAG = 8
nodes_xy = nodes
def hdist(i, j):
    (x1, y1), (x2, y2) = nodes_xy[i], nodes_xy[j]
    kx = math.cos(math.radians((y1 + y2) / 2)) * 111320
    return math.hypot((x2 - x1) * kx, (y2 - y1) * 110540)
conn = {}
for k in range(5):
    cadj = defaultdict(list)
    for a, b, _, l, _ in edges:
        if (l & 7) >= k:
            cadj[a].append(b); cadj[b].append(a)
    for (a, b), l in conn.items():
        if l >= k:
            cadj[a].append(b); cadj[b].append(a)
    comp = {}
    for s in cadj:
        if s in comp: continue
        comp[s] = s; st = [s]
        while st:
            u = st.pop()
            for v in cadj[u]:
                if v not in comp:
                    comp[v] = s; st.append(v)
    cg = defaultdict(list); C = 0.02
    for n in cadj:
        x, y = nodes_xy[n]; cg[(int(x // C), int(y // C))].append(n)
    for n in [x for x in cadj if len(cadj[x]) == 1]:
        x, y = nodes_xy[n]; gx, gy = int(x // C), int(y // C)
        best, bd = None, GAP_M
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for m in cg[(gx + dx, gy + dy)]:
                    if comp[m] != comp[n]:
                        dd = hdist(n, m)
                        if dd < bd: best, bd = m, dd
        if best is not None:
            key = (min(n, best), max(n, best))
            if conn.get(key, -1) < k: conn[key] = k
for (a, b), l in conn.items():
    edges.append((a, b, hdist(a, b), l | GAP_FLAG, []))
print('gap connectors', len(conn))

nGeom = sum(len(e[4]) for e in edges)
print('nodes', len(nodes), 'edges', len(edges), 'geom pts', nGeom)

q = lambda v: int(round(v * 1e5))
buf = bytearray(b'CVG1') + struct.pack('<III', len(nodes), len(edges), nGeom)
buf += b''.join(struct.pack('<ii', q(x), q(y)) for x, y in nodes)
gs = 0
eb, gb = bytearray(), bytearray()
for a, b, length, lvl, geom in edges:
    eb += struct.pack('<IIfB3xII', a, b, length, lvl, gs, len(geom))
    for x, y in geom:
        gb += struct.pack('<ii', q(x), q(y))
    gs += len(geom)
buf += eb + gb
with gzip.open(OUT, 'wb', compresslevel=9) as fh:
    fh.write(bytes(buf))
print('raw bytes', len(buf))

# connectivity report per vehicle class (cumulative networks)
for k in range(5):
    par = list(range(len(nodes)))
    def fd(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    used = set()
    for a, b, _, l, _ in edges:
        if (l & 7) >= k:
            par[fd(a)] = fd(b); used |= {a, b}
    from collections import Counter
    c = Counter(fd(x) for x in used)
    print('class', k, 'nodes', len(used), 'largest component %.0f%%' % (100 * max(c.values()) / len(used)))
