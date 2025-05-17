package main;

import atm.Operation;
import atm.OperationUnit;
import atm.Tries;
import org.drools.ruleunits.api.RuleUnitInstance;
import org.drools.ruleunits.api.RuleUnitProvider;

import java.util.List;

public class Main {

    public static void main(String[] args) {
        OperationUnit operationUnit = new OperationUnit();

        RuleUnitInstance<OperationUnit> instance = RuleUnitProvider.get().createRuleUnitInstance(operationUnit);

        try {

            Operation op1 = new Operation(1, true, false, true, Tries.EXPIRED, 100, 1000, 500);
            Operation op2 = new Operation(2, true, false, true, Tries.NOT_EXPIRED, 100, 1000, 500);
            Operation op3 = new Operation(3, true, false, true, Tries.NOT_EXPIRED, 1000, 100, 500);
            Operation op4 = new Operation(4, true, false, true, Tries.NOT_EXPIRED, 1000, 10000, 500);
            Operation op5 = new Operation(5, true, true, true, Tries.EXPIRED, 100, 1000, 500);
            Operation op6 = new Operation(6, false, true, true, Tries.NOT_EXPIRED, 100, 1000, 500);
            Operation op7 = new Operation(7, false, true, true, Tries.NOT_EXPIRED, 100, 1000, 500);

            operationUnit.getOperations().add(op1);
            operationUnit.getOperations().add(op2);
            operationUnit.getOperations().add(op3);
            operationUnit.getOperations().add(op4);
            operationUnit.getOperations().add(op5);
            operationUnit.getOperations().add(op6);
            operationUnit.getOperations().add(op7);

            instance.fire();
            List<Operation> authorizedOperations = instance.executeQuery("FindAuthorizedOperations").toList("$operations");
            // Repetir utilizando código JAVA

            System.out.println("AUTHORIZED: " + authorizedOperations);

        } finally {
            instance.close();
        }
    }
}
