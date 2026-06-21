# code/analyses/flu/plot_stitch.R
#   -> figures/real/flu/both_seasons/s1_nowcast_over_s2_spaghetti.pdf
# Real-time stitch (s1 nowcasts+retro top; s2 vintage spaghetti bottom).
# MechRt curves read from the girt tree; EpiNow2 from the canonical tree.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(ggtext); library(viridisLite) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R"))

col_ours <- "#2D6CA2"; col_ep <- "#5D3A1E"
windows <- list(s1 = c(as.Date("2022-10-01"), as.Date("2023-03-01")),
                s2 = c(as.Date("2023-10-01"), as.Date("2024-03-01")))
shared_x_scale <- function(w) scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y",
                                           limits = w, expand = expansion(mult = c(0.005, 0.005)))
gtree <- function(s) file.path(flu_results_dir, "girt", s)       # MechRt (girt)
ctree <- function(s) file.path(flu_results_dir, s)               # EpiNow2 (canonical)

read_mechrt_weekly <- function(s) do.call(rbind, lapply(
  list.files(file.path(gtree(s), "real_time", "results"),
             pattern = "natural_linear_taper_cv_1se", full.names = TRUE), function(f) {
    r <- readRDS(f); rt <- r$rt_df; if (is.null(rt) || !nrow(rt)) return(NULL)
    we <- as.Date(r$meta$week_end_date); rt <- rt[as.Date(rt$date) <= we, ]
    data.frame(date = as.Date(rt$date), Rt_mean = rt$Rt_mean, week_end = we) }))
read_epinow2 <- function(s) do.call(rbind, lapply(
  list.files(file.path(ctree(s), "real_time", "epinow2", "results"),
             pattern = "__gp__summary", full.names = TRUE), function(f) {
    r <- readRDS(f); rt <- r$rt_df; we <- as.Date(r$meta$week_end_date)
    rt <- rt[as.Date(rt$date) <= we & rt$type != "forecast", ]
    data.frame(date = as.Date(rt$date), Rt_mean = rt$Rt_mean, week_end = we) }))

build_top <- function() {
  lv <- c("MechRt (retrospective)","MechRt (real-time)","EpiNow2 (retrospective)","EpiNow2 (real-time)")
  retro_o <- readRDS(list.files(file.path(gtree("s1"),"retrospective","results"),
                                pattern="notaper_cv_1se", full.names=TRUE)[1])$rt_df
  retro_o <- data.frame(date=as.Date(retro_o$date), Rt_mean=retro_o$Rt_mean, series="MechRt (retrospective)")
  ow <- read_mechrt_weekly("s1")
  o_seg <- transform(ow[ow$date >= ow$week_end-6L & ow$date <= ow$week_end, ], series="MechRt (real-time)")
  o_nc  <- transform(ow[ow$date == ow$week_end, c("date","Rt_mean")], series="MechRt (real-time)")
  ep_r <- readRDS(list.files(file.path(ctree("s1"),"retrospective","epinow2","results"),
                             pattern="__gp__summary", full.names=TRUE)[1])$rt_df
  ep_r <- data.frame(date=as.Date(ep_r$date), Rt_mean=ep_r$Rt_mean, series="EpiNow2 (retrospective)")
  ew <- read_epinow2("s1")
  e_seg <- transform(ew[ew$date >= ew$week_end-6L & ew$date <= ew$week_end, ], series="EpiNow2 (real-time)")
  e_nc  <- transform(ew[ew$date == ew$week_end, c("date","Rt_mean")], series="EpiNow2 (real-time)")
  retro <- rbind(retro_o, ep_r); seg <- rbind(o_seg[,c("date","Rt_mean","week_end","series")], e_seg[,c("date","Rt_mean","week_end","series")])
  nc <- rbind(o_nc, e_nc)
  for (d in list("retro","seg","nc")) NULL
  retro$series <- factor(retro$series, levels=lv); seg$series <- factor(seg$series, levels=lv); nc$series <- factor(nc$series, levels=lv)
  pal <- c("MechRt (retrospective)"=col_ours,"MechRt (real-time)"=col_ours,"EpiNow2 (retrospective)"=col_ep,"EpiNow2 (real-time)"=col_ep)
  lty <- c("MechRt (retrospective)"="solid","MechRt (real-time)"="dashed","EpiNow2 (retrospective)"="solid","EpiNow2 (real-time)"="dotdash")
  ggplot() + geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
    geom_line(data=retro, aes(date, Rt_mean, colour=series, linetype=series), linewidth=1.0) +
    geom_line(data=seg, aes(date, Rt_mean, colour=series, linetype=series, group=interaction(series, week_end)), linewidth=1.0) +
    geom_point(data=nc, aes(date, Rt_mean, colour=series), size=1.7, show.legend=FALSE) +
    scale_colour_manual(values=pal, breaks=lv, name=NULL) + scale_linetype_manual(values=lty, breaks=lv, name=NULL) +
    shared_x_scale(windows$s1) + scale_y_continuous(breaks=seq(0,5,by=0.1)) +
    labs(x=NULL, y=expression(R[t]), title="**(a)** 2022/23 season  -  nowcasts + retrospective fits") +
    guides(colour=guide_legend(ncol=4, byrow=TRUE), linetype=guide_legend(ncol=4, byrow=TRUE)) +
    theme_minimal(base_size=14) +
    theme(plot.title=element_markdown(hjust=0.5, size=18), legend.position="bottom", legend.text=element_text(size=16),
          legend.key.width=grid::unit(2.2,"lines"), axis.text=element_text(size=14), axis.title.y=element_text(size=16),
          panel.grid.minor.y=element_blank())
}
mk_vintage_panel <- function(df, title, ylim) {
  df <- df[df$date >= windows$s2[1] & df$date <= windows$s2[2], ]
  ep_pts <- df |> group_by(week_end) |> slice_max(date, n=1L) |> ungroup()
  ggplot() + geom_hline(yintercept=1, linetype="dashed", colour="grey40") +
    geom_line(data=df, aes(date, Rt_mean, group=week_end, colour=week_end), linewidth=0.55, alpha=0.8) +
    geom_point(data=ep_pts, aes(date, Rt_mean, colour=week_end), size=1.4, alpha=0.95) +
    scale_colour_viridis_c(option="viridis", trans="date", name="Vintage",
      guide=guide_colourbar(barheight=grid::unit(7,"lines"), barwidth=grid::unit(0.6,"lines"))) +
    shared_x_scale(windows$s2) + scale_y_continuous(breaks=seq(0,5,by=0.1), limits=ylim) +
    labs(x=NULL, y=expression(R[t]), title=title) + theme_minimal(base_size=14) +
    theme(plot.title=element_markdown(hjust=0.5, size=18), axis.text=element_text(size=14), axis.title.y=element_text(size=16),
          legend.title=element_text(size=16), legend.text=element_text(size=14))
}
build_bottom <- function() {
  ep <- read_epinow2("s2"); ou <- read_mechrt_weekly("s2")
  ew <- ep[ep$date>=windows$s2[1] & ep$date<=windows$s2[2],]; ow <- ou[ou$date>=windows$s2[1] & ou$date<=windows$s2[2],]
  ylim <- range(c(ew$Rt_mean, ow$Rt_mean), na.rm=TRUE)
  (mk_vintage_panel(ep, "**(b)** 2023/24 season  -  EpiNow2", ylim) +
   mk_vintage_panel(ou, "**(c)** 2023/24 season  -  MechRt", ylim)) +
    plot_layout(guides="collect") & theme(legend.position="right")
}
stitched <- (build_top() / build_bottom()) + plot_layout(heights=c(1,1)) +
  plot_annotation(title="Real-time flu Rt across two seasons:  MechRt vs EpiNow2",
                  theme=theme(plot.title=element_text(hjust=0.5, size=18, face="bold")))
out <- file.path(flu_figures_dir, "both_seasons", "s1_nowcast_over_s2_spaghetti.pdf")
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
ggsave(out, stitched, width=13, height=8.5, dpi=300); cat("wrote", out, "\n")
