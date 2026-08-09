# Allen dorsal cortex reference

Converted from the `.npy` files shipped in [churchlandlab/wfield](https://github.com/churchlandlab/wfield)
`references/`, downloaded 260803.

## Provenance

The Allen Institute supplies only the 3D annotation volume. Everything 2-D here was
derived by wfield (`wfield/allen.py`, GPL-3, Joao Couto):

1. `mcapi.download_annotation_volume('annotation/ccf_2017', 10, ...)` -- CCF 2017 at 10 um
2. 33 dorsally visible areas listed by hand in `selection_dorsal_cortex`, resolved to
   structure ids through the adult mouse structure graph
3. `allen_top_proj_from_volume` -- first non-zero voxel seen from above, per column.
   Not a maximum projection
4. `projection_outline` -- closing, fill holes, dilation, then `find_contours(0.5)`

**Bregma is a guess.** The origin `[540 570]` is commented in wfield as "a (certainly
wrong) guess of the location of bregma". wfield aligns real data with the four
landmarks in `dorsal_cortex_landmarks.json` (OB left/centre/right at y = -3.45 mm,
RSP base at y = +3.2 mm), not by trusting this origin.

## dorsal_cortex_projection.mat

| field | size | |
|---|---|---|
| `labels` | 1320 x 1140 int16 | 0 = background, 1..33 = area |
| `x`, `y` | 1140 / 1320 | mm from `reference`, 10 um steps |
| `outline` | 4585 x 2 | closed contour of the whole map, mm |
| `reference` | [540 570] | bregma voxel, [row col] |
| `acronym`, `area_name`, `allen_id`, `allen_rgb` | 33 | from `dorsal_cortex_ccf_labels.json` |
| `area_mm2` | 33 | both hemispheres summed |

`labels` is mirror-identical about the midline, so a label alone does not name a
hemisphere. Split on `x` sign, or use the per-side contours in wfield's labels JSON.

```matlab
c = load('dorsal_cortex_projection.mat');
imagesc(c.x, c.y, c.labels); axis image; clim([0 33])
colormap([1 1 1; double(c.allen_rgb)/255])
```

## dorsal_cortex.mat

The same outline filled: `mask` 1017 x 1037 uint8 with its own `x`, `y`. It covers
77.61 mm2 against the parcellation's 76.97 -- the difference is the interhemispheric
slit at the olfactory bulbs, which the label map leaves as 0 and the fill closes.
