from matplotlib import cm, colors
def Set_palette(clu_num):
    vega_10 = list(map(colors.to_hex, cm.tab10.colors))
    vega_10_scanpy = vega_10.copy()
    # vega_10_scanpy[3] = '#ad494a'
    # vega_10_scanpy[2] = '#279e68'  # green
    # vega_10_scanpy[4] = '#aa40fc'  # purple
    vega_10_scanpy[8] = '#b5bd61'  # kakhi

    vega_20 = list(map(colors.to_hex, cm.tab20.colors))

    # reorderd, some removed, some added
    vega_20_scanpy = [
        # dark without grey:
        *vega_20[0:14:2],
        *vega_20[16::2],
        # light without grey:
        *vega_20[1:15:2],
        *vega_20[17::2],
        # manual additions:
        '#ad494a',
        '#8c6d31',
    ]
    vega_20_scanpy[2] = vega_10_scanpy[2]
    vega_20_scanpy[4] = vega_10_scanpy[4]
    vega_20_scanpy[7] = vega_10_scanpy[8]  # kakhi shifted by missing grey

    # set palette
    if clu_num <= 10:
        palette = vega_10_scanpy
        cmap = 'tab10'
    else:
        palette = vega_20_scanpy
        cmap = 'tab20'
    return palette, cmap