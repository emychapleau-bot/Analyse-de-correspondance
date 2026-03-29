# -----------------------------------------------------
# Analyse de correspondance articles de esport
# V3: Intégration de la LEMMATISATION (udpipe)
# -----------------------------------------------------

# --- INSTALLATION DES PACKAGES (si nécessaire) ---
# Si vous ne les avez pas, décommentez et exécutez ces lignes :
# install.packages("udpipe")
# install.packages("tibble")

# --- CHARGEMENT DES PACKAGES ---
library(SnowballC)
library(dplyr)
library(rlang)
library(Xplortext)
library(stopwords)

# --- NOUVEAUX PACKAGES REQUIS ---
library(udpipe)  # Le moteur de lemmatisation
library(tibble)  # Pour gérer les noms de lignes (IDs)


# ------------------------------------------------
# 1. Configuration et Chargement
# ------------------------------------------------

# 1.1. SÉLECTION DU FICHIER
filename <- file.choose()

# 1.2. CONFIGURATION DES COLONNES
cat_col <- "periodep"
text_col <- "Title"
quant_col <- "Date"
unique_id_col <- 1 

# 1.3. CHARGEMENT DES DONNÉES (CORRIGÉ)
dataFr <- read.csv2(filename, header=TRUE, row.names=1, encoding = "UTF-8")
print("Aperçu des données chargées (brutes) :")
print(head(dataFr))

# 1.4. MISE EN ORDRE
dataFr[,cat_col] <- ordered(dataFr[,cat_col], levels=c("PrePandemie", "Pandemie", "PostPandemie"))


# -----------------------------------------------------
# 2. Application de la Lemmatisation (udpipe)
# -----------------------------------------------------

# 2.1. Chargement du modèle de langue française
model_file <- "french-gsd-ud-2.5-191206.udpipe"
if (!file.exists(model_file)) {
  print("Téléchargement du modèle de langue française (french-gsd)...")
  # C'est cette fonction qui causait l'erreur
  udpipe_download_model(language = "french-gsd")
}
udmodel_french <- udpipe_load_model(file = model_file)
print("Modèle de langue chargé.")


# 2.2. Préparation des données pour l'analyse
dataFr_ids <- dataFr %>%
  tibble::rownames_to_column(var = "doc_id_temp") # Conserve les IDs (clés)

# 2.3. Lancement de la lemmatisation
data_annotated <- udpipe_annotate(udmodel_french, 
                                  x = dataFr_ids[[text_col]], 
                                  doc_id = dataFr_ids$doc_id_temp)
data_lemmas_df <- as.data.frame(data_annotated)


# 2.4. Reconstruction des titres à partir des lemmes
data_lemmatized_titles <- data_lemmas_df %>%
  filter(upos != "PUNCT") %>% 
  group_by(doc_id) %>%
  summarise(lemmatized_text = paste(tolower(lemma), collapse = " ")) %>%
  ungroup()

# 2.5. Création du dataframe final lemmatisé
dataFr_lem <- dataFr_ids %>%
  left_join(data_lemmatized_titles, by = c("doc_id_temp" = "doc_id")) %>%
  mutate(
    !!sym(text_col) := ifelse(is.na(lemmatized_text), !!sym(text_col), lemmatized_text)
  ) %>%
  select(-lemmatized_text) %>%
  tibble::column_to_rownames(var = "doc_id_temp")

print("Aperçu des données lemmatisées (lemmes) :")
print(head(dataFr_lem))


# -----------------------------------------------------
# 3. Analyse de Correspondance (Sur données lemmatisées)
# -----------------------------------------------------

# 3.1. Definition des stop words
U_stoplist<- stopwords("french") 
U_stoplist <- c(U_stoplist, "a", "e") 

# 3.2. Construction du TLE
res.TD<-TextData(dataFr_lem, var.text=text_col, idiom="fr", var.agg=cat_col, Fmin=35, Dmin=20,
                 stop.word.user=U_stoplist, stop.word.tm=TRUE, graph=FALSE) 
summary(res.TD,ndoc=0)


# MATRICE DES EFFECTIFS
TableLex<-as.matrix(res.TD$DocTerm)
dm1 <- addmargins(TableLex)
rownames(dm1)[nrow(dm1)]<-"somme"
colnames(dm1)[ncol(dm1)]<-"somme"
print(dm1, zero.print = ".")

# Matrice de profils-lignes*100
dm2<-rbind(TableLex,apply(TableLex,2,sum))
dm2<-prop.table(as.matrix(dm2),margin=1) 
dm2<-cbind(dm2,apply(dm2,1,sum))
rownames(dm2)[nrow(dm2)]<-"profil moyen"
colnames(dm2)[ncol(dm2)]<-"somme"
print(round(dm2*100,2))

# Matrice de profils-colonne*100
dm3<-cbind(TableLex,apply(TableLex,1,sum))
dm3<-prop.table(as.matrix(dm3),margin=2)
dm3<-rbind(dm3,apply(dm3,2,sum))
rownames(dm3)[nrow(dm3)]<-"somme"
colnames(dm3)[ncol(dm3)]<-"profil moyen"
round(dm3*100,2)

# calcul chi-2
res.chi2<-chisq.test(TableLex)
print(res.chi2)

# Tableaux des effectifs observés et théoriques
tabobser<-res.chi2$observed
tabobser <- addmargins(tabobser)
rownames(tabobser)[nrow(tabobser)]<-"somme"
colnames(tabobser)[ncol(tabobser)]<-"somme"
print(tabobser)

tabexpect<-res.chi2$expected
tabexpect <- addmargins(tabexpect)
rownames(tabexpect)[nrow(tabexpect)]<-"somme"
colnames(tabexpect)[ncol(tabexpect)]<-"somme"
round(tabexpect,1)

# Tableau des taux d’associations
tau<-(res.chi2$observed)/res.chi2$expected
round(tau,2)

# -----------------------------------------------------
# 4. Graphiques de l'Analyse (AFC)
# -----------------------------------------------------
res.LexCA<-LexCA(res.TD, graph=FALSE)

# graphique valeur propre
plot(res.LexCA,eigen=TRUE,selDoc=NULL,selWord=NULL,title="Valeurs propres")

# Résultats pour les documents
summary(res.LexCA,nword=0)

# Représentation des catégories d’âges
plot(res.LexCA,selWord=NULL,xlim=c(-1.9,1.9),ylim=c(-1.9,1.9),cex=1.2,col.doc="black",title="Représentation des documents")
lines(res.LexCA$row$coord[1:3,1],res.LexCA$row$coord[1:3,2],lwd=1.5,col="grey20")
lines(res.LexCA$row$coord[4:6,1],res.LexCA$row$coord[4:6,2],lwd=1.5,col="grey20")

# Représentation de tous les mots
plot(res.LexCA,gtype="DocWord",selDoc=NULL,xlim=c(-1.9,1.9),ylim=c(-1.9,1.9),col.word="black",cex=1.2,
     title="Représentation des mots")

# Les deux graphique mis ensemble
plot(res.LexCA,gtype="DocWord",xlim=c(-1.9,1.9),ylim=c(-1.9,1.9),col.doc="grey50",col.word="black",cex=1.2,
     title="Représentation simultanée des documents et des mots")
lines(res.LexCA$row$coord[1:3,1],res.LexCA$row$coord[1:3,2],lwd=1,col="black")
lines(res.LexCA$row$coord[4:6,1],res.LexCA$row$coord[4:6,2],lwd=1,col="black")

# Ellipses de confiance autour des catégories d’âge
ellipseLexCA(res.LexCA,selWord=NULL,col.doc="black",
             title="Documents - Ellipses de confiance")
lines(res.LexCA$row$coord[1:3,1],res.LexCA$row$coord[1:3,2],lwd=1,col="blue")
lines(res.LexCA$row$coord[4:6,1],res.LexCA$row$coord[4:6,2],lwd=1,col="blue")

#Ellipses de confiance autour des mots
ellipseLexCA(res.LexCA,selWord="meta 0.5",selDoc=NULL,col.word="black",
             title="Quelques mots - Ellipses de confiance")

