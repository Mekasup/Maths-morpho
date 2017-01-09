#include <iostream>

// headers pour les classes de base Image et FlatSE 
#include "../Common/Image.h"
#include "../Common/FlatSE.h"
// header pour les algos principaux de morphologie (éro, dil, ouv, fer, rec,...)
#include "../Algorithms/Morphology.h"
// header pour l'algo du Watershed
#include "../Algorithms/Watershed.h"
// header pour l'algo de seeded region growing (SRG)
#include "../Algorithms/RegionGrowing.h"
// header pour la gestion des composantes connexes
#include "../Algorithms/ConnectedComponents.h"


using namespace std;
using namespace LibTIM;

int main(int argc, char *argv[])
{
        if(argc !=2)
		{
                std::cout << "Usage: " << argv[0] << " <image.pgm>\n";
		exit(1);
                }

	// Image est une classe générique paramétrée par le type des points contenus dans l'image
        Image <unsigned char> im;
  	if(Image<U8>::load(argv[1],im))
		std::cout <<"Great, image is loaded\n";
	else return 1;


        FlatSE ptitCercle;
        ptitCercle.make2DEuclidianBall(1);

        Image<U8> imNeg = Image<U8>(im.getSize());
        for(int i=0; i<im.getSizeX(); i++)
            for(int j=0; j<im.getSizeY(); j++)
                imNeg(i,j)=255-im(i,j);
        imNeg.save("Negatif.pgm");
        Image<U8> imNegClo=closing(imNeg,ptitCercle);
        imNegClo.save("Negatif Moins Bruit");
        Image<U8> imSansFissure=erosion(imNegClo,ptitCercle);
        Image<U8> onlyFiss=imNegClo-imSansFissure;
        onlyFiss.save("Only Fissure");
        Image<U8> onlyFissEro=erosion(onlyFiss,ptitCercle);
        for(int i=0; i<im.getSizeX(); i++)
                    for(int j=0; j<im.getSizeY(); j++)
                            onlyFissEro(i,j) = onlyFissEro(i,j)*2;
        onlyFissEro.save("Only Fissure Ero");
        /*

        // Restauration
        FlatSE ligne5;
        ligneligne5.clear();
        int taille_rayon = 2;
        for(int i=-taille_rayon; i<=taille_rayon; i++) {
            Point<TCoord> point(i,0);
            ligne.addPoint(point);
        }
        ligne.setNegPosOffsets();

        Image<U8> imRest=closing(im,ligne5);
        imRest.save("Restoration.pgm");

        // Mettre en evidence les ronds blancs
        FlatSE cercle7;
        cercle7.make2DEuclidianBall(7);

        Image <U8> imOuv=opening(im, cercle7);
        imOuv.save("Ouverture.pgm");
        Image <U8> imDif=im-imOuv;
        imDif.save("Difference.pgm");

        unsigned char seuil = 125;
        for(int i=0; i<imDif.getSizeX(); i++){
            for(int j=0; j<imDif.getSizeY(); j++) {
                if(imDif(i,j)<seuil)
                    imDif(i,j) = 0;
                else
                    imDif(i,j) = 255;
            }
        }

        imDif.save("Ronds Blancs.pgm");

        //Seuillage
        Image <unsigned char> im_seuil;
        Image<U8>::load(argv[1],im_seuil);
        unsigned char seuil = 120;
        for(int i=0; i<im_seuil.getSizeX(); i++){
            for(int j=0; j<im_seuil.getSizeY(); j++) {
                if(im_seuil(i,j)<seuil)
                    im_seuil(i,j) = 0;
                else
                    im_seuil(i,j) = 255;
            }
        }

        im_seuil.save("Seuillage.pgm");

        //Erosion, dilatation
        FlatSE cercle;
        cercle.make2DEuclidianBall(3);


        Image <U8> imEro=erosion(im, cercle);
        Image <U8> imDil=dilation(im, cercle);

        imEro.save("Erosion.pgm");
        imDil.save("Dillatation.pgm");

        //Ouverture et fermeture

        Image <U8> imOpen=opening(im, cercle);
        Image <U8> imClose=closing(im, cercle);

        imOpen.save("Ouverture.pgm");
        imClose.save("Fermeture.pgm");

        //Gradient morphologique

        Image <U8> imGrad=morphologicalGradient(im, cercle);
        Image <U8> imGradInt=internalMorphologicalGradient(im, cercle);
        Image <U8> imGradFer=externalMorphologicalGradient(im, cercle);

        imGrad.save("Gradient Morpho.pgm");
        imGradInt.save("Gradient Morpho Interne.pgm");
        imGradFer.save("Gradient Morpho Externe.pgm");


        */
        /*
	// FlatSE est la classe stockant un élément structurant
	FlatSE connexity;
	// initialisation de l'élément structurant 'connexity' à un 8-voisinage (l'élément ne contient pas l'origine)
	//   X X X
	//   X . X
	//   X X X
	connexity.make2DN8();
	
	FlatSE se;
	// initialisation à un 9-voisinage (l'élément contient l'origine)
	// X X X
	// X X X
	// X X X 
	se.make2DN9();

	////////////////////////////////////////
	// exemple d'utilisation du watershed //
	////////////////////////////////////////
	
	// calcul du gradient morphologique
	Image <U8> grad=morphologicalGradient(im,se);
	grad.save("gradBeforeFiltering.pgm");
	
	// filtre h-min: suppression des minima du gradient non significatifs (évite la sur-segmentation du watershed)
	
	// paramètre h déterminant la profondeur maximale des minima éliminés
	// possibilité d'obtenir une pyramide d'images pour chaque h
	int h=atoi(argv[2]);
	
	hMinFilter(grad,se,h);
	grad.save("gradAfterFiltering.pgm");
	
	// calcul des minima régionaux du gradient filtré
	Image <U8> minima=regionalMinima(grad, connexity);
	
	// labelisation des minima
	Image <TLabel> minimaLabel1=labelConnectedComponents(minima, connexity);
	Image <TLabel> minimaLabel2=minimaLabel1;
	
	// Algo 1: watershed, algo de Meyer (le résultat est stocké dans l'image Label minimaLabel1)
	watershedMeyer(grad,minimaLabel1,connexity);
	
	// Algo 2: croissance de régions (avantage: calculé directement sur l'image originale)
	seededRegionGrowing(im,minimaLabel2,connexity);
	
	// pour chaque région labélisée, calcule la moyenne de la région dans l'image originale
	// création d'une image mosaïque
	
	Image <U8> result1=computeMarkerMean(im, minimaLabel1);
	Image <U8> result2=computeMarkerMean(im, minimaLabel2);
	
	result1.save("result1.pgm");
	result2.save("result2.pgm");
        */
}
