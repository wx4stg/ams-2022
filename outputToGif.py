from natsort import natsorted
from os import listdir, getcwd, path
import imageio

if __name__ == "__main__":
    pngs = [path.join(path.join(getcwd(), "output/"), filename) for filename in natsorted(listdir(path.join(getcwd(), "output/")))]
    with imageio.get_writer('out.gif', mode='I') as writer:
        for png in pngs:
            image = imageio.imread(png)
            writer.append_data(image)