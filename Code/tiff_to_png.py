from pathlib import Path
from PIL import Image

src = Path(r"D:\GitHub Repos\HeLa_Cell_Data\Tiffs\reyesaldasoro HeLa_Cell_Data main Tiffs-Hela_ROI_26_30_3322_2960_401")
dst = Path(r"D:\GitHub Repos\HeLa_Cell_Data\PNGs/Hela_ROI_26_30_3322_2960_401")
dst.mkdir(exist_ok=True)

for f in list(src.glob("*.tif")) + list(src.glob("*.tiff")):
    Image.open(f).save(dst / f"{f.stem}.png")
    print(f.name)