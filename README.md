# Matlab Image Crop & Rescale
This extracts any four-point-image out of any given input image array and rescales it to a defined resolution.

## Syntax
```matlab
img_out = crop_image(img_in, coordinates_x, coordinates_y, resolution)
```

`img_in`: any single color image array, if you plan on using an rgb image, you have to pass the image three times through this function.<br>
`coordinates_x`: x edge coordinates of the four point shape defining the edge of the new image. <br>
`coordinates_y`: y edge coordinates of the four point shape defining the edge of the new image.<br>
`resolution`: resolution of the new image

## Example

```matlab
img = imread('path/to/image.png');

img_cut_red = crop_image(img(:,:,1), [66, 2468, 2479, 70], [90, 100, 1436, 1460], [1920, 1080]);
img_cut_green = crop_image(img(:,:,1), [66, 2468, 2479, 70], [90, 100, 1436, 1460], [1920, 1080]);
img_cut_blue = crop_image(img(:,:,1), [66, 2468, 2479, 70], [90, 100, 1436, 1460], [1920, 1080]);
```

## Edge point order
Here's a small example of the edge points
```
   [x1, y1]
      +----------_______     [x2, y2]
      /                 --------+
     |                           \
    /                             \
   +------_________                \
 [x4, y4]          -----------------+
                                 [x3, y3]
```
