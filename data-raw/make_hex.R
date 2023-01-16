#https://www.stockio.com/free-icon/healthy-icons-crutch
library(magick)
#remotes::install_github("dmi3kno/bunny")
library(bunny)
library(ggplot2)

gst <- image_read("data-raw/design-780168de-157d-47cd-8283-0d68df47d310.png")

br_col <- "#1C1D24"
br2_col <- "#313235"
#bg_col <- "#F1EBDB"
bg_col <- "#F9F7F0"
rd_col <- "#7f7f7f"
rd1_col <- "#857E61"
rd2_col <- "#AA4465"

p <- qpd::make_pgrid(300)
ex <- -log(1-p)
rex <- log(p)
lg <- log(p)-log(1-p)

ggplot()+
  geom_line(aes(p,lg), linewidth=2, color=rd2_col)+
  geom_line(aes(p,ex), linetype=1, color=br_col)+
  annotate("text", x=0, y=1, label="Q(p)==-ln(1-p)", parse=TRUE, family="Roboto Condensed",hjust=0)+
  geom_line(aes(p,rex), linetype=1, color=br_col)+
  annotate("text", x=0.99, y=-1, label="Q(p)==ln(p)", parse=TRUE, family="Roboto Condensed",hjust=1)+
  theme_gray()+
  labs(x="p", y="x")

ggsave("data-raw/logisticQF.png", width=1550, height=1550, units="px")

plt <- image_read("data-raw/logisticQF.png")

gst_hex <- image_canvas_hex(fill_color = bg_col, border_color = br_col, border_size = 1) %>%
        image_composite(plt, gravity="center", offset="-60+0") %>%
        image_composite(gst, gravity = "center", offset = "+25-100") %>%
        image_annotate("gilchrist", gravity = "center", location = "+0+520",
                       size=200, font="Aller", color = rd_col, weight = 400) %>%
        image_annotate("gilchrist", gravity = "center", location = "+5+525",
                 size=200, font="Aller", color = br_col, weight = 400) %>%
        image_composite(image_canvas_hexborder(border_color = rd_col, border_size = 13), gravity = "center")   %>%
        image_composite(image_canvas_hexborder(border_color = br_col, border_size = 8), gravity = "center")
gst_hex %>%
  image_scale("30%")


gst_hex %>%
        image_scale("1200x1200") %>%
        image_write("data-raw/gst_hex.png", density = 600)

gst_hex %>%
        image_scale("200x200") %>%
        image_write("man/figures/logo.png", density = 600)

gst_hex_gh <- gst_hex %>%
        image_scale("400x400")


gh_logo <- bunny::github %>%
        image_scale("40x40")

gst_ghcard <- image_canvas_ghcard(fill_color = bg_col) %>%
        image_composite(gst_hex_gh, gravity = "East", offset = "+100+0") %>%
        image_annotate("Gilchrist rules", gravity = "West", location = "+100+0",
                       color=br_col, size=80, font="Aller", weight = 500) %>%
        image_annotate("Warrent G. Gilchrist (1932-2015)", gravity = "East", location = "+90+230",
                       color=br_col, size=25, font="Aller", weight = 500) %>%
        image_compose(gh_logo, gravity="West", offset = "+100+70") %>%
        image_annotate("dmi3kno/gilchrist", gravity="West",
                       location="+150+70", size=38, font="Ubuntu Mono") %>%
        image_border_ghcard(bg_col) %>%
        image_scale("50%")

gst_ghcard %>%
  image_scale("80%")

gst_ghcard %>%
        image_write("data-raw/gst_ghcard.png", density = 600)
