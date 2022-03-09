
from PIL import Image
Image.MAX_IMAGE_PIXELS = 10000000000000000000000000000000000000 ###maximum pixels
#https://stackoverflow.com/questions/10657383/stitching-photos-together

img_01= Image.open(r"F:\WRF-CHEM ANALYSIS Chap 5\rainfall\total.jpg")
img_02= Image.open(r"F:\WRF-CHEM ANALYSIS Chap 5\rainfall\direct +semi-direct.jpg")
img_03= Image.open(r"F:\WRF-CHEM ANALYSIS Chap 5\rainfall\indirect.jpg")

img_01_size = img_01.size
img_02_size = img_02.size
img_03_size = img_03.size

print('img 1 size: ', img_01_size)
print('img 2 size: ', img_02_size)
print('img 3 size: ', img_03_size)

#new_im = Image.new('RGB', size=(2 * img_01_size[0], 2 * img_01_size[1]), color=(250, 250, 250))
#https://www.geeksforgeeks.org/python-pil-image-merge-method/
new_im = Image.new('RGB',size=( (39156, 13785)))

#https://stackoverflow.com/questions/28407462/how-to-paste-an-image-onto-a-larger-image-using-pillow
new_im.paste(img_01, (0, 0))
new_im.paste(img_02, (img_01_size[0], 0))
new_im.paste(img_03, (img_01_size[0]*2,0))
#new_im.paste(img_04, (img_01_size[0], img_01_size[1]))

new_im.save(r"F:\WRF-CHEM ANALYSIS Chap 5\rainfall\merged_effects.jpg", "JPEG",quality=100)
#new_im.show()