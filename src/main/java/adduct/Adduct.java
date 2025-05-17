package adduct;


import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Adduct {


    public static int extractSign(String adduct) {
        Pattern pattern = Pattern.compile("\\[([0-9]*)M[^\\]]*](\\d*)([+-])");
        Matcher matcher = pattern.matcher(adduct);
        if (matcher.find()) {
            String chargeStr = matcher.group(2);
            return matcher.group(3).equals("-") ? -1 : 1;
        }
        return 0; // default
    }

    public static int extractCharge(String adduct) {
        if (adduct == null) return 1;

        Pattern pattern = Pattern.compile("\\[([0-9]*)M[^\\]]*](\\d*)([+-])");
        Matcher matcher = pattern.matcher(adduct);

        if (matcher.find()) {
            String chargeStr = matcher.group(2);
            return chargeStr.isEmpty() ? 1 : Integer.parseInt(chargeStr);
        }


        return 1;
    }



    public static int extractMultimer(String adduct) {
        Pattern pattern = Pattern.compile("\\[([0-9]*)M[^\\]]*](\\d*)([+-])");
        Matcher matcher = pattern.matcher(adduct);
        if (matcher.find()) {
            String multimerStr = matcher.group(1);
            return multimerStr.isEmpty() ? 1 : Integer.parseInt(multimerStr);
        }
        return 1;
    }



    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the mass difference within the tolerance respecting to the
     * massToSearch
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {
        double massToSearch;



      /*  if Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect

        if Adduct is double or triple charged the formula is M =( mz - adductMass ) * charge

        if adduct is a dimer the formula is M =  (mz - adductMass) / numberOfMultimer
*/
        if (mz != null || adduct != null) {
            int charge = extractCharge(adduct);
            int multimer = extractMultimer(adduct);
            Double adductMass = 0.0;

            if (AdductList.MAPMZPOSITIVEADDUCTS.containsKey(adduct)) {
                adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
            } else if (AdductList.MAPMZNEGATIVEADDUCTS.containsKey(adduct)) {
                adductMass = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
            }

            massToSearch = (mz * charge + adductMass) / multimer;
            return massToSearch;
        }
        return null;


}
    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {
        Double massToSearch;
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).

        /*
        if Adduct is single charge the formula is m/z = M +- adductMass. Charge is 1 so it does not affect

        if Adduct is double or triple charged the formula is mz = M/charge +- adductMass

        if adduct is a dimer or multimer the formula is mz = M * numberOfMultimer +- adductMass

        return monoisotopicMass;

         */


        if (monoisotopicMass != null || adduct != null) {
            int charge = extractCharge(adduct);
            int multimer = extractMultimer(adduct);
            Double adductMass = 0.0;

            if (AdductList.MAPMZPOSITIVEADDUCTS.containsKey(adduct)) {
                adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
            } else if (AdductList.MAPMZNEGATIVEADDUCTS.containsKey(adduct)) {
                adductMass = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
            } else {
                System.out.println(adduct + " is not a valid adduct");
                return null;
            }

            return ((monoisotopicMass * multimer) - adductMass) / charge;
        }
        return null;
        }




    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * Calcula el error entre una masa experimental y una masa teórica
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000
                / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     * Calcula la tolerance que debe haber entre el pico actual this.mz con cada peak
     * @param measuredMZ    Mz measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double measuredMZ, int ppm) {
        double deltaPPM;
        deltaPPM = (Math.abs((measuredMZ * ppm) / 1000000));
        return deltaPPM;
    }
    }




