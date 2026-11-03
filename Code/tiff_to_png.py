from pathlib import Path
from PIL import Image

src = Path(r"D:\GitHub Repos\HeLa_Cell_Data\Tiffs\Hela_ROI_30_30_2051_3944_441")
dst = Path(r"D:\GitHub Repos\HeLa_Cell_Data\PNGs/Hela_ROI_30_30_2051_3944_441")
dst.mkdir(exist_ok=True)

for f in list(src.glob("*.tif")) + list(src.glob("*.tiff")):
    Image.open(f).save(dst / f"{f.stem}.png")
    print(f.name)