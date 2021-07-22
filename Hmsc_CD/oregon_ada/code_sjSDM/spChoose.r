


spM <- otuenv %>% 
  dplyr::filter(period == "S1" & SiteName %in% train.Names) %>%
  dplyr::select(SiteName, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(SiteName, OTU) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  mutate(nSites = nSites>0) %>% # species present at M1 and M2 have value 2, this reduces this to 1
  group_by(OTU) %>% ## group just by OTU
  summarise(nSites = sum(nSites, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  tidyr::pivot_wider(names_from = trap, values_from = nSites, values_fn = function(x) sum(x)>0) %>%
  #filter(M1) %>% # CHOOOSE HERE FOR SINGLE. OR SHARED TRAP SPECIES GFROUP: filter(M1 & M2)
  # tidyr::separate(col = OTU, into = c("ID", "empty", "class", "order", "family",
  #                                         "genus", "epithet", "BOLD", "BOLDID", "size"),
  #                 remove = FALSE, sep = "_") %>%
  # dplyr::filter(order == "Diptera")%>%
  dplyr::select(OTU)

spM

sum(spM$nSites == 2)
sum(spM$nSites == 1)
unique(spM$nSites)

#####

spM <- otuenv %>% 
  dplyr::filter(period == "S1" & SiteName %in% train.Names) %>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  tidyr::pivot_wider(names_from = trap, values_from = value, values_fn = function(x) sum(x)>0) %>%
  rowwise()%>%
  mutate(M1M2 = sum(M1, M2, na.rm = T)) # change to PA


spM <- otuenv %>% 
  dplyr::filter(period == "S1" & SiteName %in% train.Names) %>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  tidyr::pivot_wider(names_from = trap, values_from = value, values_fn = function(x) sum(x)>0) %>%
  mutate(M1M2 = rowSums(select(., M1, M2), na.rm = T)>0) %>% # change to PA and if present at 1 OR 2 sites, is TRUE
  group_by(OTU) %>%
  summarise(nSites = sum(M1M2, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  
  #filter(M1) %>% # CHOOOSE HERE FOR SINGLE. OR SHARED TRAP SPECIES GFROUP: filter(M1 & M2)
  # tidyr::separate(col = OTU, into = c("ID", "empty", "class", "order", "family",
  #                                         "genus", "epithet", "BOLD", "BOLDID", "size"),
  #                 remove = FALSE, sep = "_") %>%
  # dplyr::filter(order == "Diptera")%>%
  dplyr::select(OTU)

spM
nrow(spM)