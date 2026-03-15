import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

dirname = sys.argv[1]
step = int(sys.argv[2])


def load_binary_snapshot(dirname, step):
    directory = Path(dirname)
    paths = sorted(directory.glob(f"snap*-{step:05d}.bin"))
    if not paths:
        single = directory / f"snap{step:05d}.bin"
        if not single.exists():
            raise FileNotFoundError(single)
        paths = [single]

    blocks = []
    for path in paths:
        with path.open("rb") as fp:
            time = np.fromfile(fp, np.float64, 1).item()
            nx = np.fromfile(fp, np.int32, 1).item()
            ny = np.fromfile(fp, np.int32, 1).item()
            nhyd = np.fromfile(fp, np.int32, 1).item()
            nbc = np.fromfile(fp, np.int32, 1).item()
            x = np.fromfile(fp, np.float64, nx)
            y = np.fromfile(fp, np.float64, ny)
            q = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
            b = np.fromfile(fp, np.float32, nx * ny * nbc).reshape(ny, nx, nbc)
        blocks.append((time, x, y, np.concatenate((q, b), axis=2)))

    x_starts = sorted({block[1][0] for block in blocks})
    y_starts = sorted({block[2][0] for block in blocks})
    x_offsets = {}
    y_offsets = {}
    x_parts = []
    y_parts = []
    offset = 0
    for x0 in x_starts:
        x = next(block[1] for block in blocks if block[1][0] == x0)
        x_offsets[x0] = offset
        x_parts.append(x)
        offset += x.size
    offset = 0
    for y0 in y_starts:
        y = next(block[2] for block in blocks if block[2][0] == y0)
        y_offsets[y0] = offset
        y_parts.append(y)
        offset += y.size

    x = np.concatenate(x_parts)
    y = np.concatenate(y_parts)
    fields = np.empty((y.size, x.size, blocks[0][3].shape[2]), dtype=blocks[0][3].dtype)
    time = blocks[0][0]
    for block_time, block_x, block_y, block_fields in blocks:
        if not np.isclose(block_time, time):
            raise ValueError(f"inconsistent time at step {step}")
        ix = x_offsets[block_x[0]]
        iy = y_offsets[block_y[0]]
        fields[iy:iy + block_fields.shape[0], ix:ix + block_fields.shape[1], :] = block_fields
    return x, y, time, fields

fig, ax = plt.subplots()
ax.set_xlabel("x axis")
ax.set_ylabel("y axis")

xv, yv, time, fields = load_binary_snapshot(dirname, step)
den = fields[:, :, 0]
bc = fields[:, :, 5:]
nx = xv.size
ny = yv.size

Az = np.zeros([ny,nx])

dy = yv[1] - yv[0]

for j in range(0,ny-1):
    for i in range(0,1):
        Az[j+1,i] = Az[j,i] + 0.5*( bc[j+1,i,0] + bc[j,i,0] )*dy

for j in range(0,ny):
    for i in range(0,nx-1):
        Az[j,i+1] = Az[j,i] - 0.5*( bc[j,i+1,1] + bc[j,i,1] )*dy

yy, xx = np.meshgrid(yv,xv,indexing="ij")

ax.set_xlim(xv[0], xv[-1])
ax.set_ylim(yv[0], yv[-1])
ax.text(
    0.5 * (xv[0] + xv[-1]),
    yv[-1] + 0.1 * (yv[-1] - yv[0]),
    r"$\mathrm{time}=%.2f$" % (time),
    horizontalalignment="center",
)
im = ax.imshow(
    den,
    extent=(xv[0], xv[-1], yv[0], yv[-1]),
    origin="lower",
    vmin=0,
    vmax=3,
)
ax.contour(xx, yy, Az, linestyles="solid", levels=20, colors="white")
fig.colorbar(im, ax=ax, orientation="vertical")

outputfile = dirname + "/snap%05d.pdf"%(step)

print("making plot file", outputfile)
plt.savefig(outputfile,bbox_inches="tight", pat_inches=1.0,dpi=1000)
plt.show()
