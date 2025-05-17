package atm;


import dss.Measurement;
import dss.MeasurementUnit;
import dss.RuleTest;
import org.drools.ruleunits.api.RuleUnitProvider;
import org.drools.ruleunits.api.RuleUnitInstance;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Set;

import static java.util.stream.Collectors.toList;
import static org.junit.Assert.*;


public class OperationUnitTest {


    static final Logger LOG = LoggerFactory.getLogger(OperationUnitTest.class);

    @Test
    public void testAuthorizedOperations() {
        LOG.info("Creating RuleUnit");
        OperationUnit operationUnit = new OperationUnit();

        RuleUnitInstance<OperationUnit> instance = RuleUnitProvider.get().createRuleUnitInstance(operationUnit);

        try {
            LOG.info("Insert data");
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

            LOG.info("Run query. Rules are also fired");
            instance.fire();
            List<Operation> authorizedOperations = instance.executeQuery("FindAuthorizedOperations").toList("$operations");
            assertEquals(1, authorizedOperations.size());
            assertEquals(op3, authorizedOperations.get(0));

        } finally {
            instance.close();
        }
    }

    @Test
    public void testNonVerifiedOperationRule() {
        LOG.info("Creating RuleUnit");
        OperationUnit operationUnit = new OperationUnit();

        RuleUnitInstance<OperationUnit> instance = RuleUnitProvider.get().createRuleUnitInstance(operationUnit);
        try {
            LOG.info("Insert data");
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

            LOG.info("Run query. Rules are also fired");
            instance.fire();
            Set<Operation> operaitionsNotVerified = operationUnit.getOperationsUnauthorizedNotVerified();
            assertEquals(2, operaitionsNotVerified.size());
            assertTrue(operaitionsNotVerified.contains(op6));
            assertTrue(operaitionsNotVerified.contains(op7));

            System.out.println(operationUnit.getTestString());

        } finally {
            instance.close();
        }
    }

}