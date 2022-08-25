install.packages("hexSticker")

library(hexSticker)
setwd("C:/Users/bernardo.brandao/OneDrive - NINA/presentations/Apr_2022_GIS_group_NINA")

imgurl <- "oneimpact_logo.jpg"
sticker(imgurl, package="oneimpact",
        p_size=20, p_color = "black", p_x = 1, p_y = 1.5,
        s_x=.95, s_y=.95, s_width=.85,
        h_fill = "white",
        filename="oneimpact_hex_test.png")
