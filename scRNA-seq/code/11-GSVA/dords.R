# ================================================
# Script: prepare_biomart_data.R (with Batch Processing)
# ...
# ================================================

library(biomaRt)

# ... [前面连接 BioMart 的代码保持不变] ...
# ... (假设 human_mart 和 mouse_mart 已经成功连接) ...

# 4. 获取所有小鼠基因
print("Fetching all mouse gene symbols...")
all_mouse_genes_df <- getBM(attributes = c("mgi_symbol"), mart = mouse_mart)
all_mouse_genes_vector <- all_mouse_genes_df$mgi_symbol
# 移除空值或NA
all_mouse_genes_vector <- all_mouse_genes_vector[!is.na(all_mouse_genes_vector) & all_mouse_genes_vector != ""]
print(paste("Found", length(all_mouse_genes_vector), "mouse genes to process."))


# 5. 分批次获取小鼠-人类同源基因映射
# ===================== BATCH PROCESSING LOGIC START =====================

# 定义每个批次的大小
batch_size <- 1000 
# 将基因列表分割成多个小块（chunks）
gene_chunks <- split(all_mouse_genes_vector, 
                     ceiling(seq_along(all_mouse_genes_vector) / batch_size))

print(paste("Splitting genes into", length(gene_chunks), "batches of up to", batch_size, "genes each."))

# 创建一个空的数据框来存储所有批次的结果
all_results <- data.frame()

# 使用 for 循环遍历每一个基因块
for (i in 1:length(gene_chunks)) {
  
  # 获取当前批次的基因
  current_chunk <- gene_chunks[[i]]
  
  cat(paste0("--- Processing batch ", i, " of ", length(gene_chunks), " (", length(current_chunk), " genes) ---\n"))
  
  # 使用 tryCatch 来处理单次查询可能发生的错误，使其不中断整个循环
  result_chunk <- tryCatch({
    getLDS(
      attributes = c("mgi_symbol", "ensembl_gene_id"),
      filters = "mgi_symbol",
      values = current_chunk, # <-- 关键：只查询当前批次的基因
      mart = mouse_mart,
      attributesL = c("hgnc_symbol", "ensembl_gene_id"),
      martL = human_mart,
      uniqueRows = TRUE
    )
  }, error = function(e) {
    # 如果当前批次失败，打印错误信息并返回NULL，然后继续下一个批次
    cat(paste0("--- ERROR in batch ", i, ": ", e$message, " ---\n"))
    return(NULL)
  })
  
  # 如果当前批次成功获取到结果，则将其合并到总结果中
  if (!is.null(result_chunk)) {
    all_results <- rbind(all_results, result_chunk)
    cat(paste0("--- Batch ", i, " successful. Total results so far: ", nrow(all_results), " ---\n"))
  }
  
  # 在两次请求之间暂停1秒，做一个“有礼貌的”用户，避免给服务器造成太大压力
  Sys.sleep(1)
}

# ===================== BATCH PROCESSING LOGIC END =====================


# 6. 检查并处理最终结果
if (nrow(all_results) == 0) {
  stop("FATAL ERROR: No results were returned from BioMart after all batches. Please check your query or try again later.")
}

print(paste("Successfully downloaded a map with", nrow(all_results), "total entries from all batches."))

# 7. 对数据进行预处理和清洗（代码不变）
print("Cleaning and filtering the final gene map...")
mouse_to_human_map_clean <- all_results[!is.na(all_results$HGNC.symbol) & all_results$HGNC.symbol != "", ]
mouse_to_human_map_unique <- mouse_to_human_map_clean[!duplicated(mouse_to_human_map_clean$MGI.symbol), ]
mouse_to_human_map_final <- mouse_to_human_map_unique[!duplicated(mouse_to_human_map_unique$HGNC.symbol), ]

print(paste("Final clean map contains", nrow(mouse_to_human_map_final), "unique mouse-human pairs."))

# 8. 将处理好的数据框保存到本地（代码不变）
output_dir <- "data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
file_path <- file.path(output_dir, "mouse_to_human_gene_map.rds")
saveRDS(mouse_to_human_map_final, file = file_path)

print(paste("Gene map successfully saved to:", file_path))