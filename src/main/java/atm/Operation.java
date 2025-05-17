package atm;

/**
 * @author alberto.gildelafuent
 */
public class Operation {

    private final int id_operation; // id of the operation
    private boolean verified;
    private boolean date_expired;
    private boolean correctPin; // pin of the card. Correct or not correct.
    private Tries tries; // Exceded or not exceded
    private double balance;
    private double price;
    private double limit;
    private boolean authorized; // By default to false. It will only be modified for that operations permitted by the rules.

    /**
     * @param id_operation
     * @param verified
     * @param date_expired
     * @param pin
     * @param tries
     * @param balance
     * @param price
     * @param limit
     */
    public Operation(int id_operation, boolean verified, boolean date_expired, boolean pin,
                     Tries tries, double balance, double price, double limit) {
        this.id_operation = id_operation;
        this.verified = verified;
        this.date_expired = date_expired;
        this.correctPin = pin;
        this.tries = tries;
        this.balance = balance;
        this.price = price;
        this.limit = limit;
        this.authorized = false;
    }

    public boolean isVerified() {
        return verified;
    }

    public void setVerified(boolean verified) {
        this.verified = verified;
    }

    public boolean isDate_expired() {
        return date_expired;
    }

    public void setDate_expired(boolean date_expired) {
        this.date_expired = date_expired;
    }

    public boolean isCorrectPin() {
        return correctPin;
    }

    public void setCorrectPin(boolean correctPin) {
        this.correctPin = correctPin;
    }

    public Tries getTries() {
        return tries;
    }

    public void setTries(Tries tries) {
        this.tries = tries;
    }

    public double getBalance() {
        return balance;
    }

    public void setBalance(double balance) {
        this.balance = balance;
    }

    public double getPrice() {
        return price;
    }

    public void setPrice(double price) {
        this.price = price;
    }

    public double getLimit() {
        return limit;
    }

    public void setLimit(double limit) {
        this.limit = limit;
    }

    public boolean isAuthorized() {
        return authorized;
    }

    public void setAuthorized(boolean authorized) {
        this.authorized = authorized;
    }


    @Override
    public int hashCode() {
        int hash = 7;
        hash = 89 * hash + this.id_operation;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Operation other = (Operation) obj;
        if (this.id_operation != other.id_operation) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "Operation{" + "id_operation=" + id_operation + ", verified=" + verified + ", date_expired=" + date_expired + ", pin="
                + correctPin + ", tries=" + tries + ", balance=" + balance + ", price=" + price + ", limit=" + limit + ", authorized=" + authorized + '}';
    }

}
