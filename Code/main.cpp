#include <iostream>
#include <math.h>

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

#define PI 3.14159265

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

        //Seuillage par hysterese
        FlatSE cercleHystere;
        cercleHystere.make2DEuclidianBall(5);

        int seuil_large=200;
        int seuil_etroit=151;

        Image<U8> imSeuilLarge=im;
        for(int i=0; i<im.getSizeX(); i++)
            for(int j=0; j<im.getSizeY(); j++)
                if (imSeuilLarge(i,j) < seuil_large)
                    imSeuilLarge(i,j) = 0;
        imSeuilLarge.save("Seuil Large.pgm");
        Image<U8> imSeuilEtroit=im;
        for(int i=0; i<im.getSizeX(); i++)
            for(int j=0; j<im.getSizeY(); j++)
                if (imSeuilEtroit(i,j) < seuil_etroit)
                    imSeuilEtroit(i,j) = 0;
        imSeuilEtroit.save("Seuil Etroit.pgm");

        geodesicReconstructionByDilation(imSeuilLarge,imSeuilEtroit,cercleHystere);
        imSeuilLarge.save("Seuil hystère.pgm");


        // Correction d'illumination
        FlatSE cercleIllumi;
        cercleIllumi.make2DEuclidianBall(10);

        Image <U8> riceEro=erosion(im,cercleIllumi);
        riceEro.save("Erosion.pgm");

        Image <U8> rice2=im-riceEro;
        rice2.save("Difference.pgm");

        double moyenne=0;
        int nbpixel=0;
        for(int i=0;i<riceEro.getSizeX();i++)
            for(int j=0;j<riceEro.getSizeY();j++)
            {
                nbpixel++;
                moyenne+=riceEro(i,j);
            }

        moyenne=moyenne/nbpixel;

        for(int i=0;i<rice2.getSizeX();i++)
            for(int j=0;j<rice2.getSizeY();j++)
            {
                if (rice2(i,j)<50)
                    rice2(i,j)=moyenne;
                else
                    rice2(i,j)=rice2(i,j)+moyenne;
            }
        rice2.save("RajoutFond.pgm");

        //Covariogramme angulaire
        FlatSE cercleAngulaire;
        cercleAngulaire.make2DEuclidianBall(3);
        Image <U8> imTraitement= externalMorphologicalGradient(im,cercleAngulaire);
        unsigned char seuilAngulaire = 40;
        for(int i=0; i<imTraitement.getSizeX(); i++){
            for(int j=0; j<imTraitement.getSizeY(); j++) {
                if(imTraitement(i,j)<seuilAngulaire)
                    imTraitement(i,j) = 0;
                else
                    imTraitement(i,j) = 255;
            }
        }
        imTraitement.save("TraitementSeuil.pgm");

        Image <U8> im1;
        Image <U8> im2;
        int taille_rayon_angle = 10;
        for (int ro=0; ro<180; ro++){
            FlatSE grande_ligne;
            grande_ligne.clear();
            double angle = ro;
            for(int i=-taille_rayon_angle; i<=taille_rayon_angle; i++) {
                Point<TCoord> point((int)(i*cos(angle*PI/180)),(int)(i*sin(angle*PI/180)));
                grande_ligne.addPoint(point);
            }
            grande_ligne.setNegPosOffsets();
            Image <U8> imOpenAngle=opening(imTraitement, grande_ligne);
            char nom[100];
            sprintf(nom, "3-5 Covariogramme angulaire/Seuil%d.pgm", ro);
            imOpenAngle.save(nom);
            if (ro == 61)
                im1=imOpenAngle;
            if (ro == 151)
                im2=imOpenAngle;
            int volume=0;
            for(int ii=0; ii<im.getSizeX(); ii++)
                for(int j=0; j<im.getSizeY(); j++)
                    volume += imOpenAngle(ii,j);
            cout << volume << endl;
        }

        Image<U8> imJusteTrait=im;
        for(int i=0; i<im.getSizeX(); i++)
            for(int j=0; j<im.getSizeY(); j++)
                if (im1(i,j)+im2(i,j) >= 255) {
                    imJusteTrait(i,j) = 255;
                }else{
                    imJusteTrait(i,j) = im1(i,j)+im2(i,j);
                }
        imJusteTrait.save("JusteTrait.pgm");

        // Granulometrie
        int taille = 10;
        for (int i=1; i<= taille; i++) {
            Image<U8> imDif=im;
            unsigned char seuil = 125;
            for(int ii=0; ii<imDif.getSizeX(); ii++){
                for(int jj=0; jj<imDif.getSizeY(); jj++) {
                    if(imDif(ii,jj)<seuil)
                        imDif(ii,jj) = 0;
                    else
                        imDif(ii,jj) = 255;
                }
            }
            FlatSE CercleVariable;
            CercleVariable.make2DEuclidianBall(i);
            Image<U8> imOuvVarible=opening(imDif,CercleVariable);
            char nom[100];
            sprintf(nom, "3-4 Granulometrie/Seuil%d.pgm", i);
            imOuvVarible.save(nom);
            int volume=0;
            for(int i=0; i<im.getSizeX(); i++)
                for(int j=0; j<im.getSizeY(); j++)
                    volume += imOuvVarible(i,j);
            cout << nom << volume << endl;
        }

        // Fisures
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

        unsigned char seuilRestauration = 125;
        for(int i=0; i<imDif.getSizeX(); i++){
            for(int j=0; j<imDif.getSizeY(); j++) {
                if(imDif(i,j)<seuilRestauration)
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
}
