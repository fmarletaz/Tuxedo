# ==== Packages (install once if needed) ====
# install.packages(c("ape","readr","dplyr","stringr","plotrix","svglite"))

library(ape)
library(readr)
library(dplyr)
library(stringr)
library(plotrix)  # floating.pie
library(svglite)  # for SVG export
library(tidyr)

setwd("/Users/lucas/Documents/PhD/Mespilia/Mespilia/Cafe")

# ---------- CONFIG ----------
tree_file <- "SpeciesTree_rooted_time_cal.txt"          # your ultrametric Newick in a .txt is fine
csv_file  <- "cafe2_output/total_sig_exp_contr_summary.csv"   # your table exported from Excel/CSV
pdf_out   <- "tree_node_pies.pdf"
svg_out   <- "tree_node_pies.svg"

# Colorblind-friendly palette
pie_cols <- c(
  Significant_E = "#b1494a",
  Significant_C = "#2369bd",
  No_Change     = "#bfc8d8"
)

# ---------- 1) Read tree ----------
tr <- read.tree(tree_file)
stopifnot(inherits(tr, "phylo"))
tr <- ladderize(tr, right = TRUE)   # nicer

# ---------- 2) Read CSV robustly & normalise headers ----------
raw <- read_csv(csv_file, show_col_types = FALSE)
names(raw) <- names(raw) |>
  str_replace("^\ufeff","") |>     # strip BOM
  str_trim()

# Try to find columns flexibly
find_col <- function(cands) {
  ix <- which(tolower(names(raw)) %in% tolower(cands))
  if (length(ix)) names(raw)[ix[1]] else NA_character_
}

col_node <- find_col(c("Node","node","label","Label"))
col_E    <- find_col(c("Significant_E","Significant_Expansions","Expansions","E"))
col_C    <- find_col(c("Significant_C","Significant_Contractions","Contractions","C"))
col_NC   <- find_col(c("No_Change","No Change","Nochange","Unchanged"))

if (any(is.na(c(col_node, col_E, col_C, col_NC)))) {
  stop("Could not detect required columns. Headers found: ",
       paste(names(raw), collapse=", "))
}

dat <- raw |>
  transmute(
    Node_raw = .data[[col_node]],
    Significant_E = suppressWarnings(as.numeric(.data[[col_E]])),
    Significant_C = suppressWarnings(as.numeric(.data[[col_C]])),
    No_Change     = suppressWarnings(as.numeric(.data[[col_NC]]))
  )

# Extract numeric node ids from "Holleu<1>" or "<10>"; if not present, try as integer
node_num <- suppressWarnings(as.integer(str_extract(dat$Node_raw, "(?<=<)\\d+(?=>)")))
if (all(is.na(node_num))) node_num <- suppressWarnings(as.integer(dat$Node_raw))
dat$node <- node_num

# Clean rows
dat <- dat |>
  mutate(across(c(Significant_E, Significant_C, No_Change), ~replace_na(., 0))) |>
  filter(!is.na(node)) |>
  filter(Significant_E + Significant_C + No_Change > 0)

# ---------- 3) Plot tree (base graphics) ----------
# Open devices (vector)
# --- precompute width & offsets ---
tree_width   <- max(ape::node.depth.edgelength(tr))      # total x span in branch-length units
base_r       <- 0.035 * tree_width                         # your existing pie radius
label_offset <- 2 * base_r                              # push tip labels a bit past the pies
extra_right  <- 0.06 * tree_width                         # some margin on the right for labels

# --- PDF device (Quartz on macOS) ---
pdf(pdf_out, width = 7, height = 7, family = "Helvetica", useDingbats = FALSE)
par(pty = "s", mar = c(6.5, 1.2, 1.2, 1.8))   
plot(tr,
     type = "phylogram",
     show.tip.label = TRUE,
     label.offset = label_offset,
     no.margin = TRUE,
     cex = 0.7,
     x.lim = c(0, tree_width + label_offset + extra_right))
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
radii <- rep(base_r, nrow(dat))
for (i in seq_len(nrow(dat))) {
  n  <- dat$node[i]
  x0 <- lp$xx[n]; y0 <- lp$yy[n]
  vals <- c(dat$Significant_E[i], dat$Significant_C[i], dat$No_Change[i])
  floating.pie(x0, y0, vals, radius = radii[i],
               col = c("#b1494a","#2369bd","#bfc8d8"),
               startpos = pi/2)
}


legend("topleft", inset = 0.04, bty = "n", cex = 0.8,
       fill = c("#b1494a","#2369bd","#bfc8d8"),
       legend = c("Significant expansions","Significant contractions","No change"))

# ===== Add raw counts under each species name (tips only) =====

# Settings you can tweak
tip_text_cex <- 0.7                                # size of the numbers
tip_gap_y    <- 0.7 * par("cxy")[2]                # vertical gap (in character height units)
num_cols     <- unname(pie_cols[c("Significant_E","Significant_C","No_Change")])

# helper to draw "E/C/NC" centered at x
fmt <- function(v) ifelse(abs(v - round(v)) < 1e-9, as.character(round(v)), sprintf("%.2f", v))
draw_triplet_center <- function(x, y, e, cval, nc, cex = tip_text_cex) {
  sE <- fmt(e); sC <- fmt(cval); sN <- fmt(nc)
  wE <- strwidth(sE, cex = cex); wC <- strwidth(sC, cex = cex); wN <- strwidth(sN, cex = cex)
  ws <- strwidth("/", cex = cex)
  total_w <- wE + ws + wC + ws + wN
  xstart  <- x - total_w/2
  op <- par(xpd = NA)
  text(xstart + wE/2,                y, sE, cex = cex, col = num_cols[1])
  text(xstart + wE + ws/2,           y, "/", cex = cex, col = "black")
  text(xstart + wE + ws + wC/2,      y, sC, cex = cex, col = num_cols[2])
  text(xstart + wE + ws + wC + ws/2, y, "/", cex = cex, col = "black")
  text(xstart + wE + ws + wC + ws + wN/2, y, sN, cex = cex, col = num_cols[3])
  par(op)
}

# coordinates from the tree you just plotted
lp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
Ntip <- ape::Ntip(tr)

# rows in your CSV that correspond to tips (node IDs 1..Ntip)
tip_rows <- dat %>% dplyr::filter(node <= Ntip)

# x position where tip labels start = tip's x + label_offset
x_lab <- lp$xx[tip_rows$node] + label_offset
y_lab <- lp$yy[tip_rows$node]

# center under the printed species name:
# compute tip-label width at your plotting cex (you used cex=0.7 in plot())
lab_strings <- tr$tip.label[tip_rows$node]
lab_w       <- strwidth(lab_strings, cex = 0.7)
x_center    <- x_lab + lab_w/2
y_numbers   <- y_lab - tip_gap_y

# draw the colored E/C/NC counts under each species name
for (i in seq_len(nrow(tip_rows))) {
  e  <- tip_rows$Significant_E[i]
  cc <- tip_rows$Significant_C[i]
  nc <- tip_rows$No_Change[i]
  draw_triplet_center(x_center[i], y_numbers[i], e, cc, nc, cex = tip_text_cex)
}


dev.off()
# Grab node coordinates from last plot
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
# lp$xx, lp$yy are numeric vectors indexed by node number

svglite(svg_out, width = 7, height = 7)
par(mar = c(3.5, 1, 1, 1))

plot(tr,
     type = "phylogram",
     show.tip.label = TRUE,
     label.offset = label_offset,
     no.margin = TRUE,
     cex = 0.7,
     x.lim = c(0, tree_width + label_offset + extra_right))

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
radii <- rep(base_r, nrow(dat))
for (i in seq_len(nrow(dat))) {
  n  <- dat$node[i]
  x0 <- lp$xx[n]; y0 <- lp$yy[n]
  vals <- c(dat$Significant_E[i], dat$Significant_C[i], dat$No_Change[i])
  floating.pie(x0, y0, vals, radius = radii[i],
               col = c("#b1494a","#2369bd","#bfc8d8"),
               startpos = pi/2)
}

axisPhylo(side = 1, backward = TRUE, las = 1, line = 1)

legend("topleft", inset = 0.04, bty = "n", cex = 0.8,
       fill = c("#b1494a","#2369bd","#bfc8d8"),
       legend = c("Significant expansions","Significant contractions","No change"))
dev.off()

