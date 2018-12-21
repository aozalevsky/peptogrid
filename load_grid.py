from chempy.brick import Brick
from pymol import cmd
import h5py


def load_grid(fname, prefix='aa', ramp='jet'):

    cmd.volume_ramp_new('ramp_jet', [
        0.00, 0.00, 0.00, 1.00, 0.00,
        25.00, 0.06, 0.93, 1.00, 0.05,
        50.00, 0.49, 0.51, 1.00, 0.15,
        75.00, 1.00, 1.00, 0.00, 0.23,
        90.00, 0.80, 0.18, 1.00, 0.48,
        100.00, 0.98, 0.01, 1.00, 1.00,
    ])

    cmd.volume_ramp_new('ramp_delta', [
        -1.0, 1.00, 0.00, 0.00, 0.00,
        -0.4,  1.00, 1.00, 0.00, 1.0,
        -0.2, 1.00, 0.00, 0.00, 0.00,
        0.2, 0.00, 0.00, 1.00, 0.00,
        0.4,  0.00, 1.00, 1.00, 1.0,
        1.0, 0.00, 0.00, 1.00, 0.00])

    F = h5py.File(fname, 'r')
    grid = F['step'][()]
    origin = F['origin'][()]

    NUCS = set(F.keys())
    protected = set(['origin', 'step', 'atypes'])
    NUCS -= protected

    for i in NUCS:
        data = F[i][()]
        b = Brick.from_numpy(data, grid, origin)
        bname = prefix + '_' + i
        cmd.load_brick(b, bname)
        volname = bname + '_volume'

        cmd.volume(volname, bname)

        cmd.volume_color(volname, 'ramp_jet')

# The extend command makes this runnable as a command, from PyMOL.
cmd.extend("load_grid", load_grid)
