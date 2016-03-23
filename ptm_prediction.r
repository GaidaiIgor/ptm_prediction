# ---- filename
filename = "met6_1.5.csv"
# ---- data_load
data = read.csv(filename)
# ---- data_preparation
seed = 1
filename_index = 1
formula = status ~ . -filename -type
status_ordering = c("unmodified", "modified")
# select_columns = c("filename", "status", "type", "sas", "distance", "rel_pos", "chain_size")
data = prepare_data(data)
while (TRUE)
{
  train_test = get_random_partition(data, data[, 1])
  train = train_test[[1]]
  test = train_test[[2]]
  if (nrow(train) > nrow(test))
  {
    break
  }
}
