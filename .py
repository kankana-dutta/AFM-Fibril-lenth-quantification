import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from google.colab import files
from skimage.filters import frangi
from skimage.morphology import remove_small_objects, remove_small_holes, skeletonize
from skimage.measure import label, regionprops
from scipy.ndimage import convolve

# =========================
# USER SETTINGS
# =========================
PIXEL_UM = 3.3 / 519.0   # µm per pixel

# Crop ROI to remove scalebar/colorbar (edit if needed)
Y0, Y1 = 60, 430
X0, X1 = 40, 520

# Ridge filter widths (in pixels). Tune once per dataset.
SIGMAS = np.linspace(1.0, 3.0, 7)

# Threshold on ridge-score (percentile). Higher = fewer false positives.
RIDGE_PCTL = 93.0   # try 98.5–99.5

# Morphology cleanup
MIN_OBJ_SIZE = 60
MIN_HOLE_SIZE = 60

# Keep only elongated components in the mask (before skeleton)
ECC_MIN = 0.83     # 0.88–0.97
AREA_MIN = 40       # pixels

# Length cutoff (match ImageJ)
MIN_LENGTH_UM = 0.02

# =========================
# Helper: diagonal-corrected skeleton length
# =========================
def skeleton_geodesic_length_px(component_mask: np.ndarray) -> float:
    sk = component_mask.astype(np.uint8)
    if sk.sum() < 2:
        return float(sk.sum())
    ys, xs = np.where(sk > 0)
    sset = set(zip(ys, xs))
    ortho_edges = 0
    diag_edges = 0
    forward = [(1, 0), (0, 1), (1, 1), (1, -1)]
    for y, x in sset:
        for dy, dx in forward:
            if (y + dy, x + dx) in sset:
                if dy == 0 or dx == 0:
                    ortho_edges += 1
                else:
                    diag_edges += 1
    return ortho_edges + diag_edges * np.sqrt(2.0)

# =========================
# Step 1: Upload AFM image
# =========================
uploaded = files.upload()
fname = next(iter(uploaded.keys()))
img = cv2.imread(fname, cv2.IMREAD_GRAYSCALE)
if img is None:
    raise ValueError("Could not read image.")

roi = img[Y0:Y1, X0:X1].astype(np.float32)

# =========================
# Step 2: Flatten / destripe (row-wise + col-wise background removal)
# =========================
# Row background
row_med = np.median(roi, axis=1, keepdims=True)
row_bg  = cv2.GaussianBlur(row_med, (1, 101), 0)   # smooth along y
flat = roi - row_bg

# Column background (optional but helps sometimes)
col_med = np.median(flat, axis=0, keepdims=True)
col_bg  = cv2.GaussianBlur(col_med, (101, 1), 0)   # smooth along x
flat = flat - col_bg

# Normalize to 0..1
flat = flat - flat.min()
flat = flat / (flat.max() + 1e-8)

# Mild denoise
flat_dn = cv2.GaussianBlur(flat, (3,3), 0)

# =========================
# Step 3: Ridge enhancement (Frangi)
# =========================
ridge = frangi(flat_dn, sigmas=SIGMAS, black_ridges=False)
thr = np.percentile(ridge, RIDGE_PCTL)
mask = ridge > thr

# cleanup
mask = remove_small_objects(mask, min_size=MIN_OBJ_SIZE)
mask = remove_small_holes(mask, area_threshold=MIN_HOLE_SIZE)

# =========================
# Step 4: Remove blobs; keep elongated regions
# =========================
lab = label(mask)
keep = np.zeros_like(mask, dtype=bool)

for r in regionprops(lab):
    if r.area < AREA_MIN:
        continue
    if r.eccentricity < ECC_MIN:
        continue
    keep[lab == r.label] = True

mask = keep

# =========================
# Step 5: Skeleton + topology
# =========================
skel = skeletonize(mask)

nbr = convolve(skel.astype(np.uint8), np.ones((3,3), dtype=np.uint8), mode="constant")
endpoints = skel & (nbr == 2)
branchpts = skel & (nbr >= 4)

# =========================
# Step 6: Length measurement
# =========================
lab_s = label(skel, connectivity=2)
regs = regionprops(lab_s)

rows = []
for r in regs:
    comp = (lab_s == r.label)
    Lpx = skeleton_geodesic_length_px(comp)
    Lum = Lpx * PIXEL_UM
    if Lum < MIN_LENGTH_UM:
        continue
    rows.append({"id": r.label, "length_um": float(Lum), "length_px": float(Lpx)})

df = pd.DataFrame(rows).sort_values("length_um", ascending=False).reset_index(drop=True)

# =========================
# Step 7: Report + plots
# =========================
print("\n=== Automatic AFM Fibril Report ===")
print(f"ROI: y[{Y0}:{Y1}], x[{X0}:{X1}]")
print(f"Scale: {PIXEL_UM:.6f} µm/px  (={1/PIXEL_UM:.2f} px/µm)")
print(f"Kept fibrils: {len(df)} (length > {MIN_LENGTH_UM:.2f} µm)")
if len(df):
    print(f"Mean length:   {df.length_um.mean():.3f} µm")
    print(f"Median length: {df.length_um.median():.3f} µm")
    print(f"Min / Max:     {df.length_um.min():.3f} / {df.length_um.max():.3f} µm")

fig, axs = plt.subplots(2, 3, figsize=(18, 10))
axs[0,0].imshow(roi, cmap="gray"); axs[0,0].set_title("ROI (raw)"); axs[0,0].axis("off")
axs[0,1].imshow(flat_dn, cmap="gray"); axs[0,1].set_title("Flattened + denoised"); axs[0,1].axis("off")
axs[0,2].imshow(ridge, cmap="magma"); axs[0,2].set_title(f"Frangi ridge score"); axs[0,2].axis("off")
axs[1,0].imshow(mask, cmap="gray"); axs[1,0].set_title("Mask"); axs[1,0].axis("off")
axs[1,1].imshow(skel, cmap="gray"); axs[1,1].set_title("Skeleton"); axs[1,1].axis("off")
axs[1,2].imshow(skel, cmap="gray")
by, bx = np.where(branchpts); ey, ex = np.where(endpoints)
axs[1,2].scatter(bx, by, s=10, c="red", label="branch")
axs[1,2].scatter(ex, ey, s=12, c="lime", label="end")
axs[1,2].set_title("Topology"); axs[1,2].axis("off")
axs[1,2].legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.show()

if len(df):
    plt.figure(figsize=(6,4))
    plt.hist(df["length_um"], bins=20)
    plt.xlabel("Fibril length (µm)")
    plt.ylabel("Count")
    plt.title("Length distribution")
    plt.tight_layout()
    plt.show()

# save outputs
df.to_csv("fibril_lengths_auto.csv", index=False)
cv2.imwrite("mask_auto.png", (mask.astype(np.uint8)*255))
cv2.imwrite("skeleton_auto.png", (skel.astype(np.uint8)*255))
print("Saved: fibril_lengths_auto.csv, mask_auto.png, skeleton_auto.png")

files.download("fibril_lengths_auto.csv")
files.download("mask_auto.png")
files.download("skeleton_auto.png")
