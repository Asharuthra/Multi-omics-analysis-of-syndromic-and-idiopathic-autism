library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# ----------------------------
# Read GO Sheets (Idiopathic)
# ----------------------------
bp <- read_excel("Supplementary file 3.xlsx", sheet = "Idiopathic_BP")
cc <- read_excel("Supplementary file 3.xlsx", sheet = "Idiopathic_CC")
mf <- read_excel("Supplementary file 3.xlsx", sheet = "Idiopathic_MF")

# ----------------------------
# Clean and Prepare Data
# ----------------------------
prepare_go <- function(data, ontology_name) {
  
  data %>%
    arrange(`Adjusted P-value`) %>%
    slice_head(n = 10) %>%
    mutate(
      # Remove GO ID
      Term = gsub("GO:[0-9]+", "", Term),
      
      # Remove anything inside brackets
      Term = gsub("\\(.*?\\)", "", Term),
      
      Term = str_trim(Term),
      Term = str_wrap(Term, width = 35),
      
      logP = -log10(`Adjusted P-value`),
      Ontology = ontology_name
    )
}

bp_top <- prepare_go(bp, "Biological Process")
cc_top <- prepare_go(cc, "Cellular Component")
mf_top <- prepare_go(mf, "Molecular Function")

go_all <- bind_rows(bp_top, cc_top, mf_top)

# ----------------------------
# Dark Professional Colors
# ----------------------------
custom_colors <- c(
  "Biological Process" = "#1F3B73",
  "Cellular Component" = "#2E6F57",
  "Molecular Function" = "#6A2C70"
)

# ----------------------------
# Plot (Horizontal Layout)
# ----------------------------
p <- ggplot(go_all,
            aes(x = reorder(Term, logP),
                y = logP,
                fill = Ontology)) +
  
  geom_bar(stat = "identity", width = 0.75) +
  
  coord_flip() +
  
  facet_wrap(~Ontology,
             scales = "free_y",
             nrow = 1) +
  
  scale_fill_manual(values = custom_colors) +
  
  theme_classic() +
  
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.margin = margin(15, 20, 15, 15)
  ) +
  
  labs(
    x = "GO Term",
    y = expression(-log[10]("Adjusted P-value"))
  )

print(p)

# ----------------------------
# Save High Resolution
# ----------------------------
ggsave("GO_barplot_Idiopathic_Autism.tiff",
       plot = p,
       width = 18,
       height = 7,
       dpi = 600,
       compression = "lzw")

