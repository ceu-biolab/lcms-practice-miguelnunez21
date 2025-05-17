package atm;

import java.util.HashSet;
import java.util.Set;

import dss.Measurement;
import org.drools.ruleunits.api.DataSource;
import org.drools.ruleunits.api.DataStore;
import org.drools.ruleunits.api.RuleUnitData;


public class OperationUnit implements RuleUnitData {

    private final DataStore<Operation> operations;

    private String testString;
    private final Set<Operation> operationsUnauthorizedNotVerified;

    private final Set<Operation> operationsUnauthorizedExceededLimit;


    public OperationUnit() {
        this(DataSource.createStore());
    }

    public OperationUnit(DataStore<Operation> operations) {
        this.operations = operations;
        this.operationsUnauthorizedNotVerified = new HashSet<>();
        this.operationsUnauthorizedExceededLimit = new HashSet<>();

    }

    public DataStore<Operation> getOperations() {
        return operations;
    }

    public Set<Operation> getOperationsUnauthorizedNotVerified() {
        return this.operationsUnauthorizedNotVerified;
    }

    public Set<Operation> getOperationsUnauthorizedExceededLimit() {
        return this.operationsUnauthorizedExceededLimit;
    }

    public String getTestString() {
        return testString;
    }
    public void setTestString(String testString) {
        this.testString = testString;
    }


}
