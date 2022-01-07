import matplotlib.pyplot as plt

res = []
grid_x, grid_y = 0, 0
delta_x, delta_y = 0, 0

with open('res.txt') as f:
    s = f.readline().rstrip().split()
    grid_x = int(s[0])
    grid_y = int(s[1])
    while True:
        s = f.readline()
        if not s: break
        strs = s.rstrip().split()
        vals = []
        for s in strs: vals.append(float(s))
        res.append(vals)
    
    delta_x = 1 / grid_x
    delta_y = 1 / grid_y
    
x = [delta_x * (i+1) for i in range(0, grid_x)] * grid_y
y = []
cur = delta_y
for i in range(grid_x):
    for j in range(grid_y): y.append(cur)
    cur += delta_y


#描画
fig = plt.figure(figsize = (8, 8))
ax = fig.add_subplot(111, projection='3d')
mappable = ax.scatter(y, x, res, s = 1, c=res)
fig.colorbar(mappable, shrink=0.4)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.savefig('res.png')
