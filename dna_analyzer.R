
library(ggplot2)

# DNA Sequence
dna <- readLines("sample_dna.txt")
dna <- paste(dna, collapse = "")
dna <- toupper(dna)
if (nchar(dna) == 0) {
  stop("DNA sequence is empty!")
}

if (grepl("[^ATGC]", dna)) {
  stop("Invalid DNA sequence! Only A, T, G, C allowed.")
}

print("DNA sequence is valid!")
dna_length <- nchar(dna)

cat("DNA Sequence Length:", dna_length, "\n")
gc_count <- lengths(regmatches(dna, gregexpr("[GC]", dna)))
at_count <- lengths(regmatches(dna, gregexpr("[AT]", dna)))

gc_content <- (gc_count / dna_length) * 100
at_content <- (at_count / dna_length) * 100

cat("GC Content:", round(gc_content, 2), "%\n")
cat("AT Content:", round(at_content, 2), "%\n")

# Split sequence
bases <- strsplit(dna, "")[[1]]

# Count nucleotides
counts <- data.frame(table(bases))
colnames(counts) <- c("Nucleotide", "Count")

# Calculate GC and AT content
gc <- sum(bases %in% c("G", "C"))
at <- sum(bases %in% c("A", "T"))

gc_percent <- (gc / length(bases)) * 100
at_percent <- (at / length(bases)) * 100

# Display results
print(counts)
cat("DNA Length:", length(bases), "\n")
cat("GC Content:", gc_percent, "%\n")
cat("AT Content:", at_percent, "%\n")

# Plot 1: Nucleotide frequency
p1 <- ggplot(counts, aes(x=Nucleotide, y=Count,
                         fill=Nucleotide)) +
  geom_col() +
  labs(title="DNA Nucleotide Frequency",
       x="Nucleotide", y="Count") +
  theme_minimal()

print(p1)

# Plot 2: GC vs AT content
content <- data.frame(
  Type=c("GC Content", "AT Content"),
  Percentage=c(gc_percent, at_percent)
)

p2 <- ggplot(content, aes(x="", y=Percentage,
                          fill=Type)) +
  geom_col(width=1) +
  coord_polar("y") +
  labs(title="GC vs AT Content") +
  theme_void()

print(p2)

# Save plots
ggsave("nucleotide_frequency.png", p1)
ggsave("gc_at_content.png", p2)
A_count <- sum(strsplit(dna, "")[[1]] == "A")
T_count <- sum(strsplit(dna, "")[[1]] == "T")
G_count <- sum(strsplit(dna, "")[[1]] == "G")
C_count <- sum(strsplit(dna, "")[[1]] == "C")

cat("Adenine (A):", A_count, "\n")
cat("Thymine (T):", T_count, "\n")
cat("Guanine (G):", G_count, "\n")
cat("Cytosine (C):", C_count, "\n")
