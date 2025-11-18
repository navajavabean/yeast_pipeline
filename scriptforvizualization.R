# --- 1. Load Libraries ---
library(readxl)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# --- 2. Load Data (Corrected Skip) ---
file_path <- "/Users/nmeadi00/Documents/BIOL_4315_PROJECT/eggnog_files/out.emapper.annotations.xlsx"

# If header is on line 3, we skip the first 2 lines
df <- read_excel(file_path, skip = 2)

# --- 3. Clean and Separate ---
# Force rename the first column to 'query' to be safe (sometimes it loads as '#query')
colnames(df)[1] <- "query"

# Check if GOs column exists now
if(!"GOs" %in% colnames(df)) {
  stop("Error: Still can't find 'GOs' column. Please run colnames(df) to see what the headers look like.")
}

# Process the data
df_long <- df %>%
  select(query, GOs) %>%       # Select specific columns
  filter(GOs != "-") %>%       # Remove genes with no GO terms
  separate_rows(GOs, sep = ",") # Split multiple terms (comma separated)

# --- 4. Map IDs to Terms ---
unique_gos <- unique(df_long$GOs)

go_map <- AnnotationDbi::select(GO.db,
                                keys = unique_gos,
                                columns = c("TERM", "ONTOLOGY"),
                                keytype = "GOID")

df_annotated <- left_join(df_long, go_map, by = c("GOs" = "GOID"))

# --- 5. Prepare for Plotting ---
df_annotated <- df_annotated %>%
  mutate(ONTOLOGY_FULL = case_when(
    ONTOLOGY == "BP" ~ "Biological Process",
    ONTOLOGY == "MF" ~ "Molecular Function",
    ONTOLOGY == "CC" ~ "Cellular Component",
    TRUE ~ ONTOLOGY
  ))

go_counts <- df_annotated %>%
  filter(!is.na(TERM)) %>%
  group_by(ONTOLOGY_FULL, TERM) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count)) %>%
  group_by(ONTOLOGY_FULL) %>%
  slice_head(n = 10)

# --- 6. Plot ---
ggplot(go_counts, aes(x = count, y = reorder(TERM, count), size = count, color = ONTOLOGY_FULL)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 8)) +
  facet_wrap(~ONTOLOGY_FULL, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(title = "Top GO Terms by Ontology",
       x = "Gene Count",
       y = "GO Term",
       color = "Category",
       size = "Count") +
  theme(axis.text.y = element_text(size = 8))



# ======
# --- 3. The "Ido Filter" Strategy ---
# These are the specific letters your prof wants to highlight
target_cogs <- c("E", "G", "M", "P", "I", "Q", "U", "W")

# COG categories can sometimes be "E,P" (multiple). 
# We will take the FIRST letter for the primary classification to keep it simple.
clean_df <- df %>%
  select(query, COG_category) %>%
  filter(!is.na(COG_category)) %>%
  mutate(Primary_COG = substr(COG_category, 1, 1)) # Take only the first letter

# Create the "Salt vs Other" Grouping
grouped_df <- clean_df %>%
  mutate(Category = case_when(
    Primary_COG == "E" ~ "Amino Acid Transport (E)",
    Primary_COG == "G" ~ "Carb Transport (G)",
    Primary_COG == "M" ~ "Cell Wall/Membrane (M)",
    Primary_COG == "P" ~ "Ion Transport (P)",
    Primary_COG == "I" ~ "Lipid Metabolism (I)",
    Primary_COG == "Q" ~ "Secondary Metabolites (Q)",
    Primary_COG == "U" ~ "Intracellular Trafficking (U)",
    Primary_COG == "W" ~ "Extracellular Structures (W)",
    TRUE ~ "Other" # Everything else goes here
  ))

# --- 4. Calculate Percentages ---
summary_stats <- grouped_df %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(desc(Percentage))

# --- 5. The Visualization (Bar Graph) ---
ggplot(summary_stats, aes(x = reorder(Category, Percentage), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  coord_flip() + # Flip it so it's easier to read labels
  theme_minimal() +
  labs(title = "Genomic Investment in Salt Tolerance Mechanisms",
       subtitle = "Comparison of Stress-Related COG Categories vs. Rest of Genome",
       x = "",
       y = "Percentage of Genome (%)") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), hjust = -0.1) + # Add % labels on bars
  scale_fill_brewer(palette = "Set3") + # Nice colors
  theme(legend.position = "none", # Hide legend (labels are on axis)
        axis.text.y = element_text(size = 11, face = "bold"))
