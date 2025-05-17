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
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Before;


public class OperationUnitTestWithSetUp {

    static final Logger LOG = LoggerFactory.getLogger(OperationUnitTestWithSetUp.class);

    OperationUnit operationUnit;
    RuleUnitInstance<OperationUnit> instance;
    Operation op1;
    Operation op2;
    Operation op3;
    Operation op4;
    Operation op5;
    Operation op6;
    Operation op7;

    /**
     * This method will be run before EACH single test. If we want to have common code for all tests we can use @BeforeClass
     * In this case, the creations of the facts (operations) and the insertion in each instance will be common to each test, but they should be
     * run before each test, otherwise the instance will have already trigger/fire its rules so they will not be run apart from the first test if
     * the facts do not change (refractoriness).
     */
    @Before
    public void setUp() {

        LOG.info("Creating RuleUnit");
        operationUnit = new OperationUnit();

        instance = RuleUnitProvider.get().createRuleUnitInstance(operationUnit);

        LOG.info("Insert data");
        op1 = new Operation(1, true, false, true, Tries.EXPIRED, 100, 1000, 500);
        op2 = new Operation(2, true, false, true, Tries.NOT_EXPIRED, 100, 1000, 500);
        op3 = new Operation(3, true, false, true, Tries.NOT_EXPIRED, 1000, 100, 500);
        op4 = new Operation(4, true, false, true, Tries.NOT_EXPIRED, 1000, 10000, 500);
        op5 = new Operation(5, true, true, true, Tries.EXPIRED, 100, 1000, 500);
        op6 = new Operation(6, false, true, true, Tries.NOT_EXPIRED, 100, 1000, 500);
        op7 = new Operation(7, false, true, true, Tries.NOT_EXPIRED, 100, 1000, 500);

        operationUnit.getOperations().add(op1);
        operationUnit.getOperations().add(op2);
        operationUnit.getOperations().add(op3);
        operationUnit.getOperations().add(op4);
        operationUnit.getOperations().add(op5);
        operationUnit.getOperations().add(op6);
        operationUnit.getOperations().add(op7);

        // instance.close() should not be instantiated because it will be used in each test!
    }


    @Test
    public void testAuthorizedOperations() {
        try {
            LOG.info("Run query to find authorized operations. Rules are also fired");
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
        try {
            LOG.info("Run query to find non authorized operations due to non verified cards. Rules are also fired");
            instance.fire();
            Set<Operation> operaitionsNotVerified = operationUnit.getOperationsUnauthorizedNotVerified();
            assertEquals(2, operaitionsNotVerified.size());
            assertTrue(operaitionsNotVerified.contains(op6));
            assertTrue(operaitionsNotVerified.contains(op7));

        } finally {
            instance.close();
        }
    }

}