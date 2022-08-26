# install.packages("hexSticker")
library(hexSticker)

# logo
imgurl <- "data-raw/oneimpact_logo.jpg"
# create sticker
sticker(imgurl, package = "oneimpact",
        p_size=20, p_color = "black", p_x = 1, p_y = 1.55,
        s_x=1, s_y=.95, s_width=.79,
        h_fill = "white",
        filename="man/figures/oneimpact_hex_logo.png")
