# Partie I

## 0. Préparation

1. Chargez le package `mixOmics`
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("mixOmics")
library(mixOmics)
```

2. Téléchargez et importez les données (4 fichiers: `mirna.csv`, `mrna.csv`, `protein.csv`, `sample_group.csv`)
```{r}
mirna <- read.csv("mirna.csv", header = TRUE, sep = ",", row.names = 1)
mrna <- read.csv("mrna.csv", header = TRUE, sep = ",", row.names = 1)
protein <- read.csv("protein.csv", header = TRUE, sep = ",", row.names = 1)
sample_group <- read.csv("sample_group.csv", header = TRUE, sep = ",", row.names = 1)
```

**Question 1:** Combien avez-vous d'échantillons ? de variables (mRNA, protéines, miRNA) ?
```{r}
dim(mirna)
dim(mrna)
dim(protein)
```
Les dataset contiennent tous 150 echantillons.
Pour miRNA, nous avons 184 variables.
Pour mRNA, nous avons 200 variables.
Pour protein, nous avons 142 variables.


3. Le coefficient de variation est défini comme le rapport entre l'écart-type $\sigma$ et la moyenne $\mu$ : $c_v = \frac{\sigma}{\mu}$
Construisez un fonction qui calcule le coefficient de variation à partir d'un vecteur.
```{r}
library(matrixStats)
coeff_variation <- function(dataset){
   vec_sd <- colSds(as.matrix(dataset[sapply(dataset, is.numeric)]))   # vecteur ecart-type
   vec_mean <- colMeans(as.matrix(dataset[sapply(dataset, is.numeric)]))   # vecteur moyenne
   cv <- c()  # vecteur vide
   for (i in 1:length(vec_sd)) {
      cv_tmp <- (vec_sd[i])/(vec_mean[i]) 
      cv = c(cv,cv_tmp)
   }
   return(cv)
}
```

```{r}
cv_mirna <- coeff_variation(mirna) 
cv_mrna <- coeff_variation(mrna)
cv_protein <- coeff_variation(protein)
```


4. A l'aide d'un histogramme `hist()` affichez la distribution de chacun des blocs.
```{r}
par(mfrow = c(2,2))
hist(cv_mirna, main = "Distribution du dataset 'miRNA'", xlab = NA, breaks = 25) 
hist(cv_mrna, main = "Distribution du dataset 'mRNA'", xlab = NA, breaks = 25) 
hist(cv_protein, main = "Distribution du dataset 'protein'", xlab = NA, breaks = 25)
```

**Question 2:** La distribution des coefficients de variation est-elle similaire dans les 3 blocs ?
Si oui, quel type de donnée possède le plus de variabilité ?

Nous pouvons voir que les coefficients de variation ne sont pas similaires entre les 3 blocs. 
Cependant, on voit une variabilité semblable entre les jeux de données "mirna" et "mrna" avec un intervalle de valeurs entre 0 et O,4. De plus, la distribution semble similaire et centrée.
Tandis que le jeu de données des protéines, semble varier entre -50 et 50.


5. Pour chacun des blocs, filtrez les données les plus variantes : $|c_{v}| \geq 0.15$
```{r}
mirna_filtered <- mirna[abs(cv_mirna) >= 0.15]
mrna_filtered <- mrna[abs(cv_mrna) >= 0.15]
prot_filtered <- protein[abs(cv_protein) >= 0.15]
```

**Question 3:**: Combien reste-il de gènes ? de protéines ? de miRNA ?
```{r}
cat("Le nombre de variable du dataset 'mirna' :", "\n")
cat("    dataset non-filtré", length(mirna),"variables.", "\n") # non-filtré
cat("    dataset filtré :", length(mirna_filtered),"variables.","\n") # filtré
cat("Le nombre de variable du dataset 'mrna' :", "\n")
cat("    dataset non-filtré", length(mrna),"variables.", "\n") # non-filtré
cat("    dataset filtré :", length(mrna_filtered),"variables.","\n") # filtré
cat("Le nombre de variable du dataset 'protein' :", "\n")
cat("    dataset non-filtré", length(protein),"variables.", "\n") # non-filtré
cat("    dataset filtré :", length(prot_filtered),"variables.","\n") # filtré
```
Nous pouvons voir que nous avons bien filtré les dataset de 'mirna' et de 'mrna' passant de 184 à 82 variables et de 200 à 174 variables respectivement.
Seul le dataset 'protein' n'a pas été filtré.


**Question 4:** Quel est le gène le plus variant ? La protéine associé à ce gène est-elle présente dans le jeu de donnée.
```{r}
which.max(abs(cv_mrna))
```
Le gène le plus variant est PLCD4.


**Question 5:** A l'aide des bases de donnée de votre choix répondez aux questions suivantes:

 * Quel est le rôle de ce gène ? 
Ce gène code pour un membre de la classe delta des enzymes phospholipase C jouant un rôle critique dans de nombreux processus cellulaires en hydrolysant le phosphatidylinositol 4,5-bisphosphate en deux messagers intracellulaires, l’inositol 1,4,5-trisphosphate ainsi que le diacylglycérol. L’expression de ce gène peut être un marqueur pour le cancer.
 
 * Sur quel chromosome est-il localisé ? 
Ce gène est localisé sur le chromosome 2 humain.

 * Quelle est la longueur en nucléotide de sa séquence ?
 La longueur en nucléotide de sa séquence est de 26 422 bases. (GRCh37/hg19)
    
 * Quelle est la longueur en acides aminés de la protéine associée (ou des isoformes) ?
La protéine associée au gène PLCD4 est composé de 762 acides aminés et d'ub emasse moléculaire de 87 585 Da.

\newpage

# Partie II

## 1. Single-omic: l'ACP avec `mixOmics`

**Question 6:** A quoi sert l'Analyse en Composante Principale ? Expliquez brievement sont fonctionnement ?
L'Analyse en Composante Principale permet d'analyser et de visualiser un jeu de données contenant des individus décrits par plusieurs variables quantitatives. Le principe est de diminuer le nombre de variables en créant des nouvelles variables articificielles orthogonales entre-elles permettant de visualiser un maximum de variation avec un minimum d'axe.


1. Réaliser l'ACP sur les données mRNA.
```{r}
acp <- pca(mrna_filtered, scale = FALSE, center = TRUE, ncomp = 25)  # 25 est un nombre aléatoire assez grand pour voir un bon nombre de PCA
acp
```


**Question 7:** Combien de composantes retenez-vous ? Justifiez / Illustrez
```{r}
barplot(acp$explained_variance, ylab = "Explained Variance",xlab = "Principal Components")
```
Nous retenons que les 2 premières composantes car elles représentent presque 80 % de nos données. 


2. Après avoir relancer l'ACP avec le bon nombre de composante, utiliser un graphique pour représenter les variables.
```{r}
acp2 <- pca(mrna_filtered, scale = FALSE, center = TRUE, ncomp = 2)
acp2
plotVar(acp2, cex = 2)
```

**Question 8:** Quelles sont les variables qui contribuent le plus à l'axe 1 ?
Les variables qui contribuent le plus à l'axe d'après notre pca, sont les suivants car les plus éloigné de l'axe horizontale : KDM4B, ZNF552, Cforf34, CCNA2.


3. Avec un graphique, représenter les échantillons dans l'espace formé par les composantes. 
Les échantillons sont colorés en fonction du groupe. Affichez la légende et ajoutez un titre.
```{r}
col = as.factor(sample_group$Y)
plot(acp2$x, col = col, pch = 19, main = "Analyse en Composante Principale des mRNA", xlab = "PC1(33,7%)", ylab = "PC2(23,3%)")
legend("topright", legend = c("Groupe 1", "Groupe 2", "Groupe 3"), fill = c(1,2,3), cex = 0.75)
```

4. La *sparse ACP* `spca()` implémente une étape de *feature selection*. En utilisant la documentation de la fonction et/ou l'aide disponible en ligne,  utilisez la `spca()` de manière a sélectionner 10 gènes sur la première composante et 5 gènes sur la seconde composante.
```{r}
spca <- spca(mrna_filtered, center = TRUE, scale = TRUE, ncomp = 2, keepX = c(10,5)) # 10 gènes et 5 gènes.
spca
plotVar(spca, cex = 3)
```

**Question 9:** Quelles sont les gènes que vous avez sélectionnés? *(une fonction est disponible)*
```{r}
composante1 <- selectVar(spca, comp = 1)
composante2 <- selectVar(spca, comp = 2)
composante1[1]
composante2[1]
```
Les gènes que nous avons sélectionnés sur la composante 1 sont les suivants : KDM4B, ZNF552, PREX1, TTC39A, STC2, LRIG1, C4orf34, MTL5, FUT8, SLC19A2.
Les gènes que nous avons sélectionnés sur la composante 2 sont les suivants :APBB1IP, NCF4, FLI1, C1orf162, CSF1R.  


## 2. Projection on Latent Structures

1. Réalisez une PLS `pls()` avec les données mRNA et protéines en incluant 3 composantes `(ncomp = 3)`.
```{r}
pls <- pls(X = mrna_filtered, Y = prot_filtered, ncomp = 3)
plotVar(pls, cex = c(2,2))
# plotIndiv(pls)
```
bcp prot au centre
+ de gene que de miRna

**Question 10:** A quoi sert la régression PLS pour l'intégration multi-omique?

L'idée de la régression PLS (Partial Least Squares) est de créer à partir d'un tableau de n observations décrites par p variables, un ensemble de h composantes avec h < p. 
L'avantage de cet régression est de bien s'accommoder en données manquantes. La détermination du nombre de composantes à retenir est en général fondée sur un critère mettant en jeu une validation croisée.


2. Affichez un *scatter plot* des échantillons en affichant uniquement les composantes 2 et 3.
Les échantillons doivent être coloriés par groupe. Ajoutez une légende et un titre.
```{r}
plot(acp$x[,2],acp$x[,3], main = "Analyses en Composantes Principales de PC2 et PC3", xlab = "PC2", ylab = "PC3", pch = 19, col = col)
legend("topright", legend = c("Groupe 1", "Groupe 2", "Groupe 3"), fill = c(1,2,3), cex = 0.75)
```


3. Affichez un *arrow plot* en affichant uniquement les composantes 1 et 3.
Les flèches doivent être coloriés par groupe. Ajoutez une légende et un titre.
```{r}
plotArrow(pls$x[,1],pls$x[,3]) #, main = "", xlab = "PC2", ylab = "PC3", pch = 19, col = col)
plotArrow(pls, abline = TRUE)
```


4. La *sparse PLS* `spls()` implémente une étape de *feature selection*. En utilisant la documentation de la fonction et/ou l'aide disponible en ligne,  utilisez la *sPLS* de manière a sélectionner (10 gènes, 9 protéines) sur la première composante, (5 gènes, 5 protéines) sur la seconde composante et (1 gène, 1 protéine) sur la troisième composante.
```{r}
spls <- spls(X = mrna_filtered, Y = prot_filtered, ncomp = 3, keepX = c(10,5,1), keepY = c(9,5,1))
spls
```
**Question 11:** Quels sont les variables sélectionnées sur la troisième composante.

5. Affichez un *CIM plot* à partir du modèle *sPLS*.
cim_plot_spls <- cim(spls)

**Question 12:** Quels sont les gènes et les protéines les plus corrélés? Justifiez à partir de la matrice de corrélation calculée par `cim()`.

6. Toujours à partir du même modèle *sPLS*, affichez un *network plot* en affichant uniquement les les corrélations les plus forte $(\rho \pm 0.65)$.
network_plot_spls <- network(spls)

**Question 13:** Combien de clusters / sous-graphes observés vous ?
Nous observons un seul cluster

## 2. *multiblock* Projection on Latent Structures

1. Réalisez une multiblock PLS `pls()` avec les données mRNA, protéines et miRNA `(X = list(mrna, prot), Y = mirna)` en incluant 2 composantes `(ncomp = 2)`.
multiblock_pls <- pls(X = list(mrna, protein), Y = mirna, ncomp = 2)

2. Comme la `spls()`, la `block.spls()` implémente une étape de *feature selection*. En utilisant la documentation de la fonction et/ou l'aide disponible en ligne,  utilisez la fonction de manière a sélectionner (10 gènes, 9 protéines, 7 miRNA) sur la première composante et (5 gènes, 4 protéines, 3 miRNA) sur la seconde composante.
liste <- list("mrna_filtered", "prot_filtered", "mirna_filtered")
liste2 <- list(mrna_filtered, prot_filtered, mirna_filtered)


mirna_filtered2 <- mirna_filtered
prot_filtered2 <- prot_filtered
mrna_filtered2 <- mrna_filtered

mirna_filtered2[,] <- 0
prot_filtered2[,] <- 0
mrna_filtered2 <- 0

for (k in colnames(mrna_filtered)) {

	if (colnames(mrna_filtered) == colnames(prot_filtered) & colnames(mrna_filtered) == colnames(mirna_filtered)) {
		mirna_filtered2 <- cbind(mrna_filtered[,k]) 
		prot_filtered2 <- cbind(prot_filtered[,k])
		mrna_filtered2 <- cbind(mrna_filtered[,k])
	}
}
  
block.spls(X=liste,Y=spls,ncomp=3,keepX=list(list(10,9),list(5,5),list(1,1)))
block.spls(X=liste,Y=spls,ncomp=3,keepX=list(list(10,9),list(5,5),list(1,1)))


#je ne comprends pas comment  on peut travailler par block si les noms des colonnes de blocks doivent etre les memes

**Question 14:** Quels sont les variables sélectionnées sur la première composante.

## 3. Analyse supervisée : (s)PLS-DA

Le fichier `sample_groupe.csv` associe un groupe à chaque échantillon.

**Question 15:** Donnez la répartition des groupes.

sample_group <- read.delim(file = "sample_group.csv", header = T, sep = ",")

barplot(sample_group[,1], names.arg = sample_group[,2])

plsda(X=mrna_filtered,Y= vecteur)
1. Utilisez la `pls.da()` en utilisant les gènes (`X`) et le groupe (`Y`) avec 2 composantes.

pla

2. Affichez le graphe des échantillons.

**Question 16:** Comparez ce graphe avec le graphe des échantillons obtenu avec l'ACP (1.3). Quel méthode permet d'obtenir de meilleurs clusters?

## 4. Analyse supervisée : block-(s)PLS-DA

1. Réalisez une multiblock sPLS-DA `block.splsda()` avec les données mRNA, protéines, miRNA `(X = list(mrna, prot, mirna))` et le groupe en incluant 5 composantes `(ncomp = 5)`.

2. Utiliser la fonction `perf()` sur le modèle obtenu. 

**Question 17:** Quelle serait le nombre de composante minimal à inclure ?

3. Relancez le modèle avec 2 composantes et utilisez l'option `keepX` pour sélectionner 15 gènes, protéines et miRNA sur la première compoante et 10 gènes, protéines et miRNA sur la seconde composante.

4. Réalisez un *circos plot* avec le modèle obtenu en affichant les corrélations fortes $|\rho| > 0.5$. Ajoutez un titre.

---

# Partie III

## 5. Mises en situation

Dans cette section, nous allons vous présenter deux designs expérimentaux et il
vous faudra déterminer quelle est l'approche analytique à privilégier pour
répondre aux questions demandées. Il s'agit d'intégrer à la fois l'informations
sur l'analyse bioinformatique en partant des données brutes mais également de
cibler les bonnes approches multiomiques.

1. Un de vos collègue s'intéresse aux effets de l'exposition à des polluants
   sur la santé des ours polaires. Pour ce faire, il a accès à des données
   transcriptomiques provenant d'une vingtaine de trios (un mère et sa portée
   de deux enfants) ainsi qu'à diverses mesures cliniques numériques pour tous
   les échantillons.

2. Vous travaillez sur un modèle murin et vous souhaitez comprendre les impacts
   d'un traitement sur le microbiote. Vous avez accès à des données de
   séquençage de type 16S ainsi qu'à des données de métabolomiques pour des
   souris traitées et pour des souris non-traitées. Vous pouvez prendre pour
   acquis que l'analyse primaire des données de métabolomiques a déjà été
   complétées et que vous avez déjà accès aux décomptes pour chaque molécules.
