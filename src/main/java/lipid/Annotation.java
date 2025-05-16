package lipid;

import adduct.Adduct;
import adduct.AdductList;

import java.util.*;

/**
 * Representa la anotación (feature) de un lípido.
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;          // m/z del pico escogido como cabecera
    private final double intensity;   // intensidad del pico más abundante del grupo
    private final double rtMin;       // tiempo de retención (min)
    private final Set<Peak> groupedSignals; // picos agrupados (isótopos/aductos)
    private String adduct;            // aducto inferido ([M+H]+, [M+Na]+, …)
    private double score = 0.0;
    private int totalScoresApplied = 0;
    private Ionization ionization;


    // --- constructores -----------------------------------------------------

    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, Ionization ionization ) {
        this(lipid, mz, intensity, retentionTime, Collections.emptySet(), ionization);
    }

    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, Set<Peak> groupedSignals, Ionization ionization) {

        this.lipid = lipid;
        this.mz = mz;
        this.intensity = intensity;
        this.rtMin = retentionTime;

        // Se usa TreeSet para poder acceder al pico de menor m/z rápidamente
        this.groupedSignals = orderSignals(groupedSignals);
        this.ionization = ionization;
        this.adduct = detectAdduct(this.groupedSignals);

    }

    public Set<Peak> orderSignals(Set<Peak> signals) {
        Set<Peak> orderedSignalss = new TreeSet<>(Comparator.comparingDouble(Peak::getMz));
        orderedSignalss.addAll(signals);
        return orderedSignalss;
    }


    public Lipid getLipid()      { return lipid;     }
    public double getMz()        { return mz;        }
    public double getIntensity() { return intensity; }
    public double getRtMin()     { return rtMin;     }
    public String getAdduct()    { return adduct;    }
    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }

    // --- score helpers -----------------------------------------------------

    public double getScore() { return score; }

    /** Añade un delta al score—se mantiene normalizado en [0,1]. */
    public void addScore(double delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    public double getNormalizedScore() {
        return totalScoresApplied == 0 ? 0.0 : score / totalScoresApplied;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    // ----------------------------------------------------------------------


    public String detectAdduct(Set<Peak> groupedSignals) {
        final int PPM_TOLERANCE=10;
        double deltaTolerance=Adduct.calculateDeltaPPM(this.mz,PPM_TOLERANCE);



        if (this.ionization.equals(Ionization.POSITIVE)) {
            for (String adduct1 : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
                for (String adduct2 : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
                    for (Peak p1 : groupedSignals) {
                        for (Peak p2 : groupedSignals) {
                            if (p1.equals(p2)) continue;

                            Double mono1 = Adduct.getMonoisotopicMassFromMZ(p1.getMz(), adduct1);
                            Double mono2 = Adduct.getMonoisotopicMassFromMZ(p2.getMz(), adduct2);
                            if (mono1 == null || mono2 == null) continue;

                            int ppmError=Adduct.calculatePPMIncrement(mono1,mono2);
                        if (this.lipid.getMonoisotropicmass()==0){
                            if (ppmError <= PPM_TOLERANCE){
                                if (Math.abs(p1.getMz() - this.mz) <= deltaTolerance) {
                                    return adduct1;
                                } else if (Math.abs(p2.getMz() - this.mz) <= deltaTolerance) {
                                    return adduct2;
                                }
                            }
                        }else {
                                int ppmeErrorCheck1=Adduct.calculatePPMIncrement(mono1, this.lipid.getMonoisotropicmass());
                                int ppmErrorCheck2=Adduct.calculatePPMIncrement(mono2, this.lipid.getMonoisotropicmass());

                                    if (ppmeErrorCheck1 <= PPM_TOLERANCE && ppmErrorCheck2 <= PPM_TOLERANCE) {
                                        if (Adduct.calculatePPMIncrement(p1.getMz(),this.mz) <= deltaTolerance) {
                                            return adduct1;
                                        } else if (Adduct.calculatePPMIncrement(p2.getMz(),this.mz) <= deltaTolerance) {
                                            return adduct2;
                                        }
                                    }
                                }

                            }
                    }
                }
            }

        } else if (this.ionization.equals(Ionization.NEGATIVE)) {
            for (String adduct1 : AdductList.MAPMZNEGATIVEADDUCTS.keySet()) {
                for (String adduct2 : AdductList.MAPMZNEGATIVEADDUCTS.keySet()) {
                    for (Peak p1 : groupedSignals) {
                        for (Peak p2 : groupedSignals) {
                            if (p1.equals(p2)) continue;

                            Double mono1 = Adduct.getMonoisotopicMassFromMZ(p1.getMz(), adduct1);
                            Double mono2 = Adduct.getMonoisotopicMassFromMZ(p2.getMz(), adduct2);
                            if (mono1 == null || mono2 == null) continue;

                            int ppmError = Adduct.calculatePPMIncrement(mono1, mono2);
                            if (this.lipid.getMonoisotropicmass() == 0) {
                                if (ppmError <= PPM_TOLERANCE) {
                                    if (Math.abs(p1.getMz() - this.mz) <= deltaTolerance) {
                                        return adduct1;
                                    } else if (Math.abs(p2.getMz() - this.mz) <= deltaTolerance) {
                                        return adduct2;
                                    }
                                }
                            } else {
                                int ppmeErrorCheck1 = Adduct.calculatePPMIncrement(mono1, this.lipid.getMonoisotropicmass());
                                int ppmErrorCheck2 = Adduct.calculatePPMIncrement(mono2, this.lipid.getMonoisotropicmass());

                                if (ppmeErrorCheck1 <= PPM_TOLERANCE && ppmErrorCheck2 <= PPM_TOLERANCE) {
                                    if (Adduct.calculatePPMIncrement(p1.getMz(), this.mz) <= deltaTolerance) {
                                        return adduct1;
                                    } else if (Adduct.calculatePPMIncrement(p2.getMz(), this.mz) <= deltaTolerance) {
                                        return adduct2;
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
        return "unknown";
    }





    public static void printScores(Annotation... annotations) {
        for (int i = 0; i < annotations.length; i++) {
            System.out.printf("Final score of annotation%d: %.1f%n", i + 1, annotations[i].getNormalizedScore());
        }
        System.out.println("---------------------------------------------------------------");
    }



    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation a)) return false;
        return Double.compare(a.mz, mz) == 0 &&
                Double.compare(a.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, a.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format(
                "Annotation(%s, mz=%.5f, RT=%.2f, adduct=%s, intensity=%.1f, score=%.3f)",
                lipid.getName(), mz, rtMin, adduct, intensity, score
        );
    }
}
