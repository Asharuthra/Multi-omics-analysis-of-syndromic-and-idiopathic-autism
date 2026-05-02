library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# ---------------------------------
# Read KEGG Sheets (NEW FILE NAME)
# ---------------------------------
syn_kegg  <- read_excel("Kegg.xlsx", sheet = "Syndromic_KEGG")
idio_kegg <- read_excel("Kegg.xlsx", sheet = "Idiopathic_KEGG")

# ---------------------------------
# Prepare Function
# ---------------------------------
prepare_kegg <- function(data, group_name) {
  
  data %>%
    arrange(`Adjusted P-value`) %>%
    slice_head(n = 15) %>%
    mutate(
      # Extract gene count from Overlap (e.g., 8/63 → 8)
      Count = as.numeric(sub("/.*", "", Overlap)),
      
      # Clean pathway name
      Term = gsub("\\(.*?\\)", "", Term),
      Term = str_trim(Term),
      Term = str_wrap(Term, width = 40),
      
      logP = -log10(`Adjusted P-value`),
      Group = group_name
    )
}

syn_top  <- prepare_kegg(syn_kegg,  "Syndromic Autism")
idio_top <- prepare_kegg(idio_kegg, "Idiopathic Autism")

kegg_all <- bind_rows(syn_top, idio_top)

# ---------------------------------
# Bubble Plot (Separate Panels)
# ---------------------------------
p <- ggplot(kegg_all,
            aes(x = logP,
                y = reorder(Term, logP),
                size = Count,
                color = logP)) +
  
  geom_point(alpha = 0.85) +
  
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  
  scale_color_gradient(low = "#6BAED6", high = "#08306B") +
  
  theme_classic() +
  
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold")
  ) +
  
  labs(
    x = expression(-log[10]("Adjusted P-value")),
    y = "KEGG Pathway",
    size = "Gene Count",
    color = expression(-log[10]("Adj P"))
  )

print(p)

# ---------------------------------
# Save High Resolution
# ---------------------------------
ggsave("KEGG_Bubble_Separate_Panels.tiff",
       plot = p,
       width = 16,
       height = 8,
       dpi = 600,
       compression = "lzw")
