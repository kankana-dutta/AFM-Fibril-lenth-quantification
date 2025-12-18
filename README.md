# Automatic AFM Fibril Length Quantification 

A fully automated, reproducible pipeline to detect AFM fibrils and measure their lengths (µm) from a single grayscale AFM image. The workflow is designed to minimize manual interaction while producing ImageJ-like length outputs plus QC figures and intermediate masks.

---

## What this program does

Given an AFM image, the script:

1. **Loads the image and crops an ROI**  
   Removes scale bars / color bars / borders so downstream processing focuses only on the physical sample region.

2. **Flattens the AFM background and removes striping**  
   AFM images often contain line-by-line offsets and slow background variation.  
   This code estimates a **row-wise** and **column-wise** background using medians + smoothing, then subtracts them.  
   Result: fibrils become higher-contrast and more uniformly detectable across the field.

3. **Enhances fibril-like ridges (Frangi filter)**  
   Uses a multi-scale ridge detector (Frangi) to boost elongated, tube-like structures while suppressing noise and blobs.

4. **Thresholds + morphological cleanup**  
   Converts the ridge score into a binary fibril mask using a percentile threshold, then removes tiny objects and small holes.

5. **Rejects non-fibril blobs by shape**  
   Keeps only elongated connected components using **eccentricity** and **area** filtering.

6. **Skeletonizes the fibril mask**  
   Reduces each fibril to a 1-pixel-wide centerline.

7. **Measures fibril lengths with diagonal correction**  
   Computes a geodesic skeleton length using 4-neighbor + diagonal edges (diagonal edges weighted by √2), then converts pixels → µm using the provided calibration.

8. **Exports results + QC outputs**
   - `fibril_lengths_auto.csv` (sorted lengths)
   - `mask_auto.png` (binary fibril mask)
   - `skeleton_auto.png` (skeleton overlay source)
   - Summary stats + diagnostic plots (raw ROI, flattened image, ridge map, mask, skeleton, endpoints/branch points)

---

## Why these steps are necessary (AFM-specific)

AFM datasets commonly include:
- **Row/scan-line artifacts** and **slow background curvature**  
- **Spatially varying contrast**
- Noise that can fragment thin fibrils or create false positives

Flattening + ridge enhancement makes detection robust and consistent, so length measurements reflect fibril geometry rather than imaging artifacts.

---

## Outputs

- **CSV:** `fibril_lengths_auto.csv`  
  Columns: `id`, `length_um`, `length_px`
- **Images:**  
  `mask_auto.png`, `skeleton_auto.png`
- **Plots:**  
  ROI, flattened image, ridge score, mask, skeleton + topology markers, and length histogram.

---

## Key parameters to tune

- `Y0, Y1, X0, X1` : ROI crop bounds  
- `RIDGE_PCTL` : ridge-score threshold (higher → fewer detections, fewer false positives)  
- `SIGMAS` : ridge widths (match fibril thickness range)  
- `ECC_MIN`, `AREA_MIN` : removes non-fibril objects  
- `MIN_LENGTH_UM` : discard tiny fragments / noise
---

 ## Raw image  
<img width="300" height="300" alt="Screenshot 2025-12-18 180222" src="https://github.com/user-attachments/assets/2932bc88-da99-4e16-972d-333b7f9290fa" />

## Processed image
<img width="2466" height="678" alt="image" src="https://github.com/user-attachments/assets/09f12c26-47af-4bcd-b7e9-23dbdd123fdc" />
<img width="2456" height="676" alt="image" src="https://github.com/user-attachments/assets/584b2c9d-3dfa-4107-aa5f-71b5b792901d" />

## Estimated Results:
=== Automatic AFM Fibril Report ===
- ROI: y[60:430], x[40:520]
- Scale: 0.006358 µm/px  (=157.27 px/µm)
- Kept fibrils: 31 (length > 0.02 µm)
- Mean length:   0.234 µm
- Median length: 0.192 µm
- Min / Max:     0.093 / 0.543 µm
---

### Note
Compared to manual ImageJ quantification, the summary statistics produced by this script may still carry an estimated ~5% error; the pipeline is designed to reduce and standardize this error via consistent preprocessing, ridge detection, skeletonization, and length measurement.

