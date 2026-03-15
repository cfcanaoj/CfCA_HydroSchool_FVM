import sys
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

dirname = sys.argv[1]
step_s = int(sys.argv[2])
step_e = int(sys.argv[3])


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

fname_anime = "animation.mp4"

graph_list = []
for istep in range(step_s, step_e + 1):
    print("making plot ", f"{dirname}/snap*-{istep:05d}.bin")
    x, y, time, fields = load_binary_snapshot(dirname, istep)
    den = fields[:, :, 0]

    xmin = x[0]
    xmax = x[-1]
    ymin = y[0]
    ymax = y[-1]
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    pg00 = ax.text(
        0.5 * (xmin + xmax),
        ymax + 0.1 * (ymax - ymin),
        r"$\mathrm{time}=%.3f$" % (time),
        horizontalalignment="center",
    )
    im = ax.imshow(
        den,
        extent=(xmin, xmax, ymin, ymax),
        origin="lower",
        vmin=0,
        vmax=3,
        animated=True,
    )

    if istep == step_s:
        fig.colorbar(im, ax=ax, orientation="vertical")
    graph_list.append([pg00, im])


ani = animation.ArtistAnimation(fig, graph_list, interval=200)
print("making animation file", fname_anime)
ani.save(dirname + "/" + fname_anime, writer="imagemagick")
plt.show()
