data_list <- gnomesims::gnome_mx_simulation(npgsloci = 10)

gnome_power_data <- data_list$power %>%
  mutate(across(nmz:ndz, ~ round(., 0))) %>%
  mutate(across(a:x, ~ round(., 2))) %>%
  mutate(across(p1:p8, ~ round(., 3))) %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz)) %>%
  dplyr::rename(`CT` = g, `SI` = b, `CT(m1) MZDZ` = p1,
         `SI(m2) MZDZ` = p2, `CT(m3) MZDZ` = p3, `SI(m3) MZDZ` = p4,
         `CT(m1) DZ` = p5, `SI(m2) DZ` = p6, `CT(m3) DZ` = p7,
         `SI(m3) DZ` = p8)


gnome_params_data <- data_list$params %>%
  mutate(across(nmz:ndz, ~ round(., 0))) %>%
  mutate(across(a:x, ~ round(., 2))) %>%
  mutate(across(e1:e8, ~ round(., 3))) %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz)) %>%
  dplyr::rename(`CT` = g, `SI` = b, `CT(m1) MZDZ` = e1,
                `SI(m2) MZDZ` = e2, `CT(m3) MZDZ` = e3, `SI(m3) MZDZ` = e4,
                `CT(m1) DZ` = e5, `SI(m2) DZ` = e6, `CT(m3) DZ` = e7,
                `SI(m3) DZ` = e8)

save(gnome_power_data, file = "~/Documents/GitHub/gnomesims/data/gnome_power_data.RData")
save(gnome_params_data, file = "~/Documents/GitHub/gnomesims/data/gnome_params_data.RData")

library(usethis)

# Read the dataset
load("~/Documents/GitHub/gnomesims/data/gnome_params_data.RData")
load("~/Documents/GitHub/gnomesims/data/gnome_power_data.RData")

# Save the processed dataset in the `data/` directory
usethis::use_data(gnome_power_data, overwrite = TRUE)
usethis::use_data(gnome_params_data, overwrite = TRUE)
