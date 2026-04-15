import io
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow

# Pillow is required for the saving/compression routine below.
# If not installed: pip install pillow
from PIL import Image

# CRISPR-Cas9 লক্ষ্যমাত্রার সাইট, কাট এবং এইচডিআরের জন্য ভিজ্যুয়ালাইজেশন
def visualize_crispr_workflow(dna_sequence, targets, cut_position=None, repair_template=None):
    print("Input DNA Sequence:")
    print(dna_sequence)
    print("\nCRISPR Target Sites (Guide RNA + PAM):")
    for i, (position, guide, pam) in enumerate(targets, start=1):
        print(f"  Target {i}:")
        print(f"    Position: {position}")
        print(f"    Guide RNA: {guide}")
        print(f"    PAM: {pam}")
    if cut_position is not None:
        print(f"\nCut Position: {cut_position}")
    if repair_template is not None:
        print(f"\nHDR Repair Template: {repair_template}")

    # ভিজুয়ালাইজেশন শুরু
    fig, ax = plt.subplots(figsize=(len(dna_sequence) // 2, 6))  # ডিএনএ সিকোয়েন্সের দৈর্ঘ্যের উপর ভিত্তি করে প্রস্থ গতিশীলভাবে সমন্বয়

    nucleotide_colors = {
        'A': 'blue',
        'T': 'yellow',
        'G': 'green',
        'C': 'red'
    }

    # ডিএনএ চেনের ভিজ্যুয়ালাইজেশন (শীর্ষ সারি)
    for i, base in enumerate(dna_sequence):
        color = nucleotide_colors[base]
        ax.add_patch(Rectangle((i, 3), 1, 1, color=color, alpha=0.7))  # উচ্চতা ৩ এ অবস্থান

    # ড্যাশযুক্ত আয়তক্ষেত্র দিয়ে CRISPR লক্ষ্যমাত্রার সাইটগুলি (গাইড RNA + PAM) চিহ্নিত
    for position, guide, pam in targets:
        ax.add_patch(Rectangle((position, 2), len(guide), 1, fill=False, edgecolor='black', linewidth=1.5, linestyle='--'))  # গাইড RNA
        ax.add_patch(Rectangle((position + len(guide), 2), len(pam), 1, fill=False, edgecolor='orange', linewidth=1.5, linestyle='--'))  # PAM

    # সিমুলেটেড কাট ভিজুয়ালাইজেশন (যদি কাট অবস্থান দেওয়া থাকে)
    if cut_position is not None:
        ax.axvline(x=cut_position, ymin=0.5, ymax=2.5, color='red', linestyle='-', linewidth=1.5, label="Cut Position")
        ax.text(cut_position + 0.2, 2.7, 'Cut', color='red', fontsize=10, weight='bold')

    # এইচডিআর ভিজ্যুয়ালাইজেশন (যদি মেরামত টেম্পলেট দেওয়া হয়)
    if repair_template is not None:
        repair_start = cut_position
        repair_end = cut_position + len(repair_template)
        for i, base in enumerate(repair_template):
            color = nucleotide_colors[base]
            ax.add_patch(Rectangle((repair_start + i, 1), 1, 1, color=color, alpha=0.7, edgecolor='black', linewidth=1.5))
        ax.text(repair_start + len(repair_template) / 2, 0.7, 'HDR Template', fontsize=10, color='black', weight='bold', ha='center')

    # সিকোয়েন্সের জন্য লেবেল যুক্ত
    plt.text(len(dna_sequence) / 3, 4.2, "Input DNA Sequence", fontsize=12, color="black", weight="bold")
    plt.text(len(dna_sequence) / 3, 2.2, "Target DNA Sequence with PAM", fontsize=12, color="black", weight="bold")

    # বিস্তারিত তথ্যের জন্য লেজেন্ড যোগ
    legend_patches = [
        Rectangle((0, 0), 1, 1, color=color, alpha=0.7, label=base) for base, color in nucleotide_colors.items()
    ] + [
        Rectangle((0, 0), 1, 1, edgecolor='black', fill=False, linestyle='--', label="Guide RNA Target"),
        Rectangle((0, 0), 1, 1, edgecolor='orange', fill=False, linestyle='--', label="PAM Sequence"),
        Rectangle((0, 0), 1, 1, edgecolor='red', facecolor='none', label="Cut Position"),
    ]
    ax.legend(handles=legend_patches, loc='lower left', title="Notations", fontsize=8)

    # চিত্রের নান্দনিকতা সামঞ্জস্য
    ax.set_xlim(0, len(dna_sequence))
    ax.set_ylim(0, 5)
    ax.axis('off')  # একটি ক্লিনার ভিজ্যুয়ালাইজেশনের জন্য এক্সিসগুলি বন্ধ
    plt.title("CRISPR-Cas9 Workflow Visualization")

    # ----------------- Saving routine: TIFF at 1200 DPI and < 5 MB -----------------
    # Strategy summary:
    # - Render figure to an in-memory PNG at DPI=1200 to obtain the pixel representation.
    # - Convert to PIL.Image and attempt to save as TIFF with JPEG compression (tiff_jpeg).
    # - If file size > 5MB, reduce JPEG quality iteratively. If still too large, downsample pixels
    #   (while keeping DPI metadata set to 1200) and retry.
    #
    # IMPORTANT: DPI metadata in the final TIFF will be set to 1200 regardless of pixel downsampling.
    out_filename = "crispr_visualization_1200dpi.tiff"
    size_limit_bytes = 5 * 1024 * 1024
    REQUIRED_DPI = 1200

    def fig_to_pil_at_dpi(fig_obj, dpi):
        buf = io.BytesIO()
        fig_obj.savefig(buf, format="png", dpi=dpi, bbox_inches='tight', pad_inches=0)
        buf.seek(0)
        img = Image.open(buf).convert('RGB')
        buf.close()
        return img

    # Render at the requested DPI to start
    pil_img = fig_to_pil_at_dpi(fig, REQUIRED_DPI)
    orig_size = pil_img.size  # (width, height) in pixels
    # Attempt to save with JPEG-in-TIFF varying quality
    def tiff_bytes_from_pil(img, quality):
        out = io.BytesIO()
        try:
            img.save(out, format="TIFF", compression="tiff_jpeg", quality=quality, dpi=(REQUIRED_DPI, REQUIRED_DPI))
            data = out.getvalue()
            out.close()
            return data
        except Exception:
            out.close()
            return None

    saved = False
    final_bytes = None
    final_scale = 1.0
    final_quality = None

    # Try quality sweep on original pixels first
    for q in range(95, 4, -5):  # 95,90,...,5
        data = tiff_bytes_from_pil(pil_img, q)
        if data is None:
            # if tiff_jpeg not supported by Pillow build, break and try fallback compressions/downsampling
            break
        if len(data) <= size_limit_bytes:
            final_bytes = data
            final_quality = q
            saved = True
            break

    # If still not saved under limit, progressively downsample and repeat quality loop.
    if not saved:
        # downsample multipliers (reduce pixel count while keeping DPI metadata at 1200)
        scales = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15]
        for scale in scales:
            new_w = max(1, int(orig_size[0] * scale))
            new_h = max(1, int(orig_size[1] * scale))
            resized = pil_img.resize((new_w, new_h), resample=Image.LANCZOS)
            for q in range(95, 4, -5):
                data = tiff_bytes_from_pil(resized, q)
                if data is None:
                    break
                if len(data) <= size_limit_bytes:
                    final_bytes = data
                    final_quality = q
                    final_scale = scale
                    saved = True
                    break
            if saved:
                break

    # Fallbacks if tiff_jpeg unsupported or not successful:
    if not saved:
        # Try other TIFF compressions (lossless) — some may produce smaller files for vector-like plots.
        fallback_compressions = ["tiff_adobe_deflate", "tiff_deflate", "tiff_lzw"]
        for comp in fallback_compressions:
            out = io.BytesIO()
            try:
                pil_img.save(out, format="TIFF", compression=comp, dpi=(REQUIRED_DPI, REQUIRED_DPI))
                data = out.getvalue()
                out.close()
                if len(data) <= size_limit_bytes:
                    final_bytes = data
                    final_quality = None
                    final_scale = 1.0
                    saved = True
                    break
            except Exception:
                out.close()
                continue

    # Final aggressive fallback: downsample heavily and save with low JPEG quality, still embed DPI=1200
    if not saved:
        aggressive_scales = [0.12, 0.10, 0.08]  # extreme downsampling as last resort
        for scale in aggressive_scales:
            new_w = max(1, int(orig_size[0] * scale))
            new_h = max(1, int(orig_size[1] * scale))
            small = pil_img.resize((new_w, new_h), resample=Image.LANCZOS)
            data = tiff_bytes_from_pil(small, 5)  # min quality
            if data is not None and len(data) <= size_limit_bytes:
                final_bytes = data
                final_quality = 5
                final_scale = scale
                saved = True
                break

    # If still not saved under limit, save a best-effort TIFF (smallest produced) or PNG fallback
    if not saved:
        # If we produced any TIFF bytes during attempts, pick the smallest and save it.
        # Otherwise save a PNG fallback (PNG may be smaller for such plots).
        candidate_path = os.path.splitext(out_filename)[0] + "_best_effort.tiff"
        smallest = None
        # Try capturing one last tiff_jpeg at very low quality from a small downsample
        try:
            small = pil_img.resize((max(1, int(orig_size[0]*0.15)), max(1, int(orig_size[1]*0.15))),
                                   resample=Image.LANCZOS)
            data = tiff_bytes_from_pil(small, 5)
            if data is not None:
                smallest = data
        except Exception:
            smallest = None

        if smallest is not None:
            with open(candidate_path, "wb") as f:
                f.write(smallest)
            print(f"Saved best-effort TIFF (may be >5MB): {candidate_path} ({len(smallest)/(1024*1024):.2f} MB)")
            print("Note: could not achieve <5MB with available compressions without unacceptable downsampling.")
        else:
            # PNG fallback
            png_path = os.path.splitext(out_filename)[0] + ".png"
            try:
                pil_img.save(png_path, format='PNG', dpi=(REQUIRED_DPI, REQUIRED_DPI), optimize=True)
                print(f"Saved PNG fallback: {png_path} ({os.path.getsize(png_path)/(1024*1024):.2f} MB). "
                      "TIFF <5MB not achievable with available compressions.")
            except Exception as e:
                print("Final save fallback failed:", str(e))
    else:
        # Write the successful TIFF bytes to disk and report details
        with open(out_filename, "wb") as f:
            f.write(final_bytes)
        size_mb = len(final_bytes) / (1024 * 1024)
        print(f"\nSaved TIFF: {out_filename}")
        print(f"Size: {size_mb:.2f} MB (limit 5.00 MB)")
        print(f"DPI metadata: {REQUIRED_DPI}. Pixel scale factor applied: {final_scale:.3f}. JPEG quality used: {final_quality}")

    # ----------------- End saving routine -----------------
    plt.show()

dna_seq = "TCCAGGATCGGGTGGGATCTCATTGTTCAA"
targets = [(0, "GGTGGGATCTCATTGTTCAA", "AGG")]
cut_position = 17  # উদাহরণ কাটার স্থান
repair_template = "GTTTCCTGCCACAGGGTCATGCTCTTTAAGCTCTCAGAAGAAGTGAGCGAGTTGGA"  # উদাহরণ HDR টেমপ্লেট

visualize_crispr_workflow(dna_seq, targets, cut_position, repair_template)
