package lipid;

import java.util.Objects;

public class Lipid {
    private final int compoundId;
    private final String name;
    private final String formula;
    private final LipidType lipidType;
    private final int carbonCount;
    private final int doubleBondsCount;
    private double monoisotropicmass;


    /**
     * @param compoundId
     * @param name
     * @param formula
     * @param lipidType
     * @param carbonCount
     * @param doubleBondCount
     */
    public Lipid(int compoundId, String name, String formula, LipidType lipidType, int carbonCount, int doubleBondCount) {
        this.compoundId = compoundId;
        this.name = name;
        this.formula = formula;
        this.lipidType = lipidType;
        this.carbonCount = carbonCount;
        this.doubleBondsCount = doubleBondCount;
    }
    public Lipid(int compoundId, String name, String formula, LipidType lipidType, int carbonCount, int doubleBondCount, double monoisotropicmass) {
        this.compoundId = compoundId;
        this.name = name;
        this.formula = formula;
        this.lipidType = lipidType;
        this.carbonCount = carbonCount;
        this.doubleBondsCount = doubleBondCount;
        this.monoisotropicmass=monoisotropicmass;
    }



    public int getCompoundId() {
        return compoundId;
    }

    public String getName() {
        return name;
    }

    public String getFormula() {
        return formula;
    }

    public LipidType getLipidType() {
        return lipidType;
    }

    public int getCarbonCount() {
        return carbonCount;
    }

    public int getDoubleBondsCount() {
        return doubleBondsCount;
    }

    public double getMonoisotropicmass() {
        return monoisotropicmass;
    }

    public void setMonoisotropicmass(double monoisotropicmass) {
        this.monoisotropicmass = monoisotropicmass;
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof Lipid lipid)) return false;
        return compoundId == lipid.compoundId;
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(compoundId);
    }

    @Override
    public String toString() {
        return "Lipid{" +
                "compoundId=" + compoundId +
                ", name='" + name + '\'' +
                ", formula='" + formula + '\'' +
                ", lipidType='" + lipidType + '\'' +
                ", carbonCount=" + carbonCount +
                ", doubleBondCount=" + doubleBondsCount +
                '}';
    }
}
