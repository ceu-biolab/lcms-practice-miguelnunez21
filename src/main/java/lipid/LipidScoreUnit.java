package lipid;

import org.drools.ruleunits.api.DataSource;
import org.drools.ruleunits.api.DataStore;
import org.drools.ruleunits.api.RuleUnitData;


import java.util.Map;

public class LipidScoreUnit implements RuleUnitData {
    private final DataStore<Annotation> annotations;
    public LipidScoreUnit() {
        this(DataSource.createStore());
    }

    public LipidScoreUnit(DataStore<Annotation> annotations) {
        this.annotations = annotations;

    }


    public DataStore<Annotation> getAnnotations() {
        return annotations;
    }


    private static final Map<LipidType, Integer> PRIORITY = Map.of(
            LipidType.PG, 1,
            LipidType.PE, 2,
            LipidType.PI, 3,
            LipidType.PA, 4,
            LipidType.PS, 5,
            LipidType.PC, 10
    );

    public static boolean hasHigherRT(LipidType type1, LipidType type2) {
        Integer p1 = PRIORITY.get(type1);
        Integer p2 = PRIORITY.get(type2);
        if (p1 == null || p2 == null) return false;
        return p1 > p2;
    }
}


