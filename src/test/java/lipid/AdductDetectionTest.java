package lipid;


import adduct.Adduct;
import org.junit.Before;
import org.junit.Test;

import lipid.Annotation;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;



public class AdductDetectionTest {
    // !!TODO For the adduct detection both regular algorithms or drools can be used as far the tests are passed.


    @Before
    public void setup() {
        // !! TODO Empty by now,you can create common objects for all tests.
    }

    @Test
    public void shouldDetectAdductBasedOnMzDifference() {
        // [M+H]+ (700.500)  vs [M+Na]+ (~21.98 Da más)
        Peak mH  = new Peak(700.500, 100000.0);
        Peak mNa = new Peak(722.482,  80000.0);
        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", LipidType.TG, 34, 1);

        Annotation annotation = new Annotation(
                lipid, 700.49999, 80000.0, 6.5, Set.of(mH, mNa), Ionization.POSITIVE);




        assertNotNull("El aducto no debería ser null", annotation.getAdduct());
        assertEquals("[M+H]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectLossOfWaterAdduct() {
        Peak mh = new Peak(700.500, 90000.0);        // [M+H]+
        Peak mhH2O = new Peak(682.4894, 70000.0);     // [M+H–H₂O]+, ~18.0106 Da less

        Lipid lipid = new Lipid(1, "PE 36:2", "C41H78NO8P", LipidType.TG, 36, 2);
        Annotation annotation = new Annotation(lipid, mhH2O.getMz(), mhH2O.getIntensity(), 7.5d, Set.of(mh, mhH2O), Ionization.POSITIVE);





        assertNotNull("[M+H-H2O]+ should be detected", annotation.getAdduct());

        assertEquals( "Adduct inferred from lowest mz in group","[M+H-H2O]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectDoublyChargedAdduct() {
        // Assume real M = (700.500 - 1.0073) = 699.4927
        // So [M+2H]2+ = (M-(- 2.0146))/2 = 350.75365
        Peak singlyCharged = new Peak(700.500, 100000.0);  // [M+H]+
        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+


        Lipid lipid = new Lipid(3, "TG 54:3", "C57H104O6", LipidType.TG, 54, 3);
        Annotation annotation = new Annotation(lipid, doublyCharged.getMz(), doublyCharged.getIntensity(), 10d, Set.of(singlyCharged, doublyCharged), Ionization.POSITIVE);

        assertNotNull("[M+2H]2+ should be detected", annotation.getAdduct());

        assertEquals( "[M+2H]2+", annotation.getAdduct());
    }




        @Test
        public void testDoublyChargedAdduct() {
            // [M+2H]2+ → M = (mz*charge + adductMass)
            String adduct = "[M+2H]2+";
            double mz = 350.7536;  // Valor típico
            double expectedNeutralMass = (mz*2+(-1.007276*2)); // Base de la fórmula
            System.out.println(expectedNeutralMass);

            Double actualNeutralMass = Adduct.getMonoisotopicMassFromMZ(mz, adduct);

            assertNotNull("Neutral mass should not be null", actualNeutralMass);
            System.out.println("Neutral mass: " + actualNeutralMass);
            assertEquals("Calculated neutral mass should match expected", expectedNeutralMass, actualNeutralMass, 0.0001);
        }

    @Test
    public void testMZFromMonoisotopicMass_PositiveAdducts() {
        double neutralMass = 700.0;

        // [M+H]+ => m/z = (M / 1) + adductMassl
        Double mzH = Adduct.getMZFromMonoisotopicMass(neutralMass, "[M+H]+");
        System.out.println(mzH);
        assertNotNull(mzH);
        assertEquals("Check [M+H]+", (neutralMass / 1) - (-1.007276), mzH, 0.0001);

        // [M+Na]+
        Double mzNa = Adduct.getMZFromMonoisotopicMass(neutralMass, "[M+Na]+");
        System.out.println(mzNa);
        assertNotNull(mzNa);
        assertEquals("Check [M+Na]+", (neutralMass / 1) - (-22.989218), mzNa, 0.0001);

        // [M+2H]2+ => m/z = ((M / 2) + adductMass)
        Double mz2H = Adduct.getMZFromMonoisotopicMass(neutralMass, "[M+2H]2+");
        assertNotNull(mz2H);
        assertEquals("Check [M+2H]2+", ((neutralMass ) - (-1.007276*2))/2, mz2H, 0.0001);

        // [2M+H]+ => m/z = ((M * 2) / 1) + adductMass
        Double mz2M = Adduct.getMZFromMonoisotopicMass(neutralMass, "[2M+H]+");
        assertNotNull(mz2M);
        assertEquals("Check [2M+H]+", (neutralMass * 2 / 1) - (-1.007276), mz2M, 0.0001);
    }

    @Test
    public void shouldDetectNegativeAdductBasedOnMzDifference() {
        // [M‑H]− (698.485)  vs  [M+Cl]− (734.462)
        Peak mHminus = new Peak(698.485, 100000.0);
        Peak mCl     = new Peak(734.462,  80000.0);

        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", LipidType.TG, 34, 1);

        Annotation annotation = new Annotation(lipid, mHminus.getMz(), mHminus.getIntensity(), 6.5, Set.of(mHminus, mCl), Ionization.NEGATIVE);



        assertNotNull("El aducto no debería ser null", annotation.getAdduct());
        assertEquals("[M-H]-", annotation.getAdduct());
    }

    @Test
    public void shouldDetectDoublyChargedNegativeAdduct() {
        // Simulamos un lípido con masa real aproximada (M)
        // [M−H]− con m/z = 699.485
        // [M+2H]2− ≈ (M - 2.0146) / 2 = ~349.238862

        Peak singlyCharged = new Peak(699.485, 100000.0);  // [M−H]−
        Peak doublyCharged = new Peak(349.239, 85000.0);   // [M+2H]2− (hipotético, debe estar en MAPMZNEGATIVEADDUCTS)

        Lipid lipid = new Lipid(4, "PE 36:2", "C41H78NO8P", LipidType.TG, 36, 2);

        Annotation annotation = new Annotation(lipid, doublyCharged.getMz(), doublyCharged.getIntensity(), 9.2, Set.of(singlyCharged, doublyCharged), Ionization.NEGATIVE);
             System.out.println(annotation.getAdduct());

        assertNotNull("Debería detectar un aducto", annotation.getAdduct());
        assertEquals("Debería detectar [M-2H]2- como aducto base", "[M-2H]2-", annotation.getAdduct());
    }

    @Test
    public void shouldDetectDimerizationAndBaseMonomer() {
        // Monómero [M+H]+ 700.500
        double monoisotropicmass=699.492724;

        Peak monomer = new Peak(700.500, 100000.0);
        // Dímero [2M+H]+ 1400.993 (2*699.492724)-(-1.007276)=1399.993
        Peak dimer = new Peak(1399.993, 50000.0);

        Lipid lipid = new Lipid(5, "PC 34:2", "C42H80NO8P", LipidType.TG, 34, 2, monoisotropicmass);

        Annotation annotation = new Annotation(lipid, monomer.getMz(), monomer.getIntensity(), 7.2, Set.of(monomer, dimer), Ionization.POSITIVE);
        //concepto malo poniendo en annotationmz el dimero para obtener el monomero

        // Verificamos que se detecta correctamente como [M+H]+ (el más pequeño)
        assertNotNull("Debe detectar el aducto principal", annotation.getAdduct());
        assertEquals("Debe ser [M+H]+ como base", "[M+H]+", annotation.getAdduct());
    }


    //-----------------------------------------------------------------------




    @Test
    public void shouldDetectClAdduct() {
        // Simulamos una molécula con M = 700.500
        // El aducto [M+Cl]− sería: (M + 34.969402) = 735.469402
        // mapMZNegativeAdductsTMP.put("[M-H]−", 1.007276d);
        // mapMZNegativeAdductsTMP.put("[M+Cl]−", -735.469402d);

        Peak mH = new Peak(698.4854, 100000.0); // [M–H]⁻
        Peak mCl = new Peak(734.462, 80000.0);  // [M+Cl]⁻, diferencia de ~34.9694 Da      // [M-H]− (más bajo, pero menor intensidad)

        Lipid lipid = new Lipid(3, "TG 54:3", "C57H104O6", LipidType.TG, 54, 3);
        Annotation annotation = new Annotation(lipid, mH.getMz(), mH.getIntensity(), 10d, Set.of(mH, mCl), Ionization.NEGATIVE);


        assertNotNull("[M-H]- should be detected", annotation.getAdduct());
        assertEquals("[M-H]-", annotation.getAdduct());
    }

    @Test
    public void shouldDetectMostProbableAdductAmongFormateAndMH() {
        double neutralMass = 699.492724;
        // Simulamos una molécula con M = 700.500
        // [M-H]− => 700.500 - 1.007276 = 699.492724
        // [M+HCOOH-H]− => 699.492724 + 44.998201 = 744.490925

        Peak mhAdductPeak = new Peak(neutralMass - 1.007276d, 80000.0);
        Peak formateAdductPeak = new Peak(neutralMass + 44.998201 , 85000.0);

        Lipid lipid = new Lipid(5, "PE 38:4", "C43H80NO8P", LipidType.PE, 38, 4);
        Annotation annotation = new Annotation(lipid, mhAdductPeak.getMz(), mhAdductPeak.getIntensity(), 10d,  Set.of(mhAdductPeak, formateAdductPeak),Ionization.NEGATIVE);

        assertNotNull("Adduct should not be null", annotation.getAdduct());
        assertEquals("[M-H]-", annotation.getAdduct());
    }

    @Test
    public void shouldDetectBestAdductAmongThreeNegativeOptions() {
        //Masa neutra simulada de la molécula
        double neutralMass = 699.492724;

        Peak clAdductPeak = new Peak(neutralMass + 34.969402, 80000.0);
        Peak formateAdductPeak = new Peak(neutralMass + 44.998201, 90000.0);
        Peak mhH2OAdductPeak = new Peak(neutralMass - 1.0073 - 18.0106, 85000.0);
        Peak doublyChargedPeak = new Peak((neutralMass - 1.007276*2)/2, 75000.0);
        //Peak mhAdductPeak = new Peak(neutralMass - 1.0073, 80000.0);

        Lipid lipid = new Lipid(6, "TG 34:1", "C39H79N2O6P", LipidType.TG, 34, 1);
        Annotation annotation = new Annotation(lipid, clAdductPeak.getMz(), clAdductPeak.getIntensity(), 10d, Set.of(clAdductPeak, formateAdductPeak, mhH2OAdductPeak, doublyChargedPeak), Ionization.NEGATIVE);


        assertNotNull("Adduct should be inferred from available peaks", annotation.getAdduct());
        assertEquals("[M+Cl]-", annotation.getAdduct());
    }



    @Test
    public void shouldDetectAdductFromMultiplePeaks() {
        Peak p1 = new Peak(700.500, 100000.0);// [M+H]+
        Peak p2 = new Peak(722.482, 80000.0);// [M+Na]+
        Peak p3 = new Peak(350.754, 85000.0);// [M+2H]2+
        Peak p4 = new Peak(682.4894, 70000.0);// [M+H–H₂O]+
        
        Lipid lipid = new Lipid(5, "PC 36:4", "C44H80NO8P", LipidType.TG, 36, 4);


        Annotation annotation = new Annotation(lipid, p1.getMz(), p1.getIntensity(), 6.0,Set.of(p1, p2, p3, p4),Ionization.POSITIVE );
        assertEquals("[M+H]+", annotation.getAdduct());
        assertEquals("[M+H]+", annotation.getAdduct());  // Es el que coincide con annotationMz
    }


}



