Bonjour!

Il serait utile de préciser que lors de la réalisation de notre projet, nous avons remarqué que:

Selon si l’on utilise les minuties de la première liste (m1) en tant que centre de rotation, ou les minuties de la deuxième liste (m2) en tant que centre de rotation, nous retrouvons des résultats de comparaison différents.

En effet, d’après la réponse obtenue sur Piazza : https://piazza.com/class/ktijhp746sr283?cid=148

Nous avons jugé correct d’utiliser en tant que centre de rotation les minuties de la deuxième liste (m2), et donc, appliquer la transformation sur les minuties de la première liste (m1), ce qui semble être plus correct si l’on suit l’instruction donnée dans l’énoncé “Écrivez enfin une méthode applyTransformation qui applique dans l’ordre une rotation et une translation à une minutie donnée.” puisque la position des minuties m2 ne change pas, comme indiqué dans la réponse.

Nous obtenons alors les mêmes résultats de comparaison que ceux qui sont donnés par Mme. Barbara Jobstmann dans le fichier texte results_center_m2.txt. Notamment, pour la comparaison de l’image 1_1.png avec 1_6.png nous avons aussi un faux négatif.

De plus, nous avons ajouté une fonctionnalité à notre méthode match() qui affiche le nombre de superpositions entre chacune des images: Match Cnt (Première valeur).

Cependant, ce dernier s’avère être légèrement différent de celui donné dans le fichier texte.





Afin de modulariser le code, nous avons créé deux méthodes suplémentaires dans Fingerprint:
1. findBoolean (ligne 75) permettant de savoir si un pixel au coordonnées [row, col] est noir ou blanc.
    La fonction retourne true si le pixel à ces coordonnées est noir. Elle retourne false si le pixel est blanc, ou si il se trouve en dehors des bornes de l'image.

2. copy (ligne 194) permettant de faire une copie du tableau bidimensionnel (de l'image).






Explication de l'algorithme pour la méthode ConnectedPixels (ligne 238):

Tout d'abord, nous créons un tableau de booleans à deux dimensions vide, appelé pixRes (ce sera le tableau avec uniquement les pixels connectés).
Puis, nous colorions le pixel correspondant au coordonées de la minutie (pixRes[row][col]) en noir, directement (car la minutie est en effet connectée à elle-même).
Nous avons alors au moins un pixel noir, déjà présent dans l'image resultante. 
C'est grace à cet unique pixel, que nous pourrons affirmer plus tard qu'un autre pixel y est connecté.

Ensuite, nous initialisons une variable booléenne pixChanged, grace à la quelle nous saurons quand faut-il sortir de la boucle do-while.
En entrant dans la boucle, on met la valeur de cette variable à false, afin de pouvoir sortir de la boucle si aucun changement n'a été fait.

C'est alors que nous parcourons chaque pixel de l'image de départ, en vérifiant si ce pixel satisfait les conditions suivantes:
- le pixel est noir
- il y a aumoins un voisin noir du pixel à ces coordonnées dans l'image résultante
- le pixel ce trouve dans le carré donné par la distance en paramètre

Si un pixel satisfait ces trois conditions, alors on regarde dans l'image résultante, si ce pixel était blanc initialement.

Si c'est le cas, alors on l'ajoute à l'image résultante, et on change la valeur du boolean pixChanged en true, puisque l'on vient d'effectuer un changement.

Dès-que tout les pixels connectés à la minutie ont été ajoutés dans pixRes, pixChanged n'est plus évalué à true. Par conséquent, on sort de la boucle.

Enfin, on retourne le tableau des pixels connectés.

