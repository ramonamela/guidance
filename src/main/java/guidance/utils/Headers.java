package guidance.utils;

import java.util.HashMap;

import guidance.GuidanceImpl;


public class Headers {

    /**
     * Construct the X Header
     * 
     * @return
     */
    public static String constructHeaderX(String sex) {
        StringBuilder headerMixedXBuilder = new StringBuilder();
        headerMixedXBuilder.append("chr").append("\t");
        headerMixedXBuilder.append("position").append("\t");
        headerMixedXBuilder.append("rs_id_all").append("\t");
        headerMixedXBuilder.append("info_all").append("\t");
        headerMixedXBuilder.append("certainty_all").append("\t");
        headerMixedXBuilder.append("alleleA").append("\t");
        headerMixedXBuilder.append("alleleB").append("\t");
        headerMixedXBuilder.append("all_A").append("\t");
        headerMixedXBuilder.append("all_B").append("\t");
        headerMixedXBuilder.append("all_AA").append("\t");
        headerMixedXBuilder.append("all_AB").append("\t");
        headerMixedXBuilder.append("all_BB").append("\t");
        headerMixedXBuilder.append("all_NULL").append("\t");
        headerMixedXBuilder.append("all_total").append("\t");
        headerMixedXBuilder.append("all_maf").append("\t");
        headerMixedXBuilder.append("all_info").append("\t");
        headerMixedXBuilder.append("all_impute_info").append("\t");
        headerMixedXBuilder.append("cases_A").append("\t");
        headerMixedXBuilder.append("cases_B").append("\t");
        headerMixedXBuilder.append("cases_AA").append("\t");
        headerMixedXBuilder.append("cases_AB").append("\t");
        headerMixedXBuilder.append("cases_BB").append("\t");
        headerMixedXBuilder.append("cases_NULL").append("\t");
        headerMixedXBuilder.append("cases_total").append("\t");
        headerMixedXBuilder.append("cases_maf").append("\t");
        headerMixedXBuilder.append("cases_info").append("\t");
        headerMixedXBuilder.append("cases_impute_info").append("\t");
        headerMixedXBuilder.append("controls_A").append("\t");
        headerMixedXBuilder.append("controls_B").append("\t");
        headerMixedXBuilder.append("controls_AA ").append("\t");
        headerMixedXBuilder.append("controls_AB ").append("\t");
        headerMixedXBuilder.append("controls_BB ").append("\t");
        headerMixedXBuilder.append("controls_NULL").append("\t");
        headerMixedXBuilder.append("controls_total").append("\t");
        headerMixedXBuilder.append("controls_maf").append("\t");
        headerMixedXBuilder.append("controls_info").append("\t");
        headerMixedXBuilder.append("controls_impute_info").append("\t");
        if(GuidanceImpl.getSex1().equals(sex)) {
            headerMixedXBuilder.append("sex=1_A ").append("\t");
            headerMixedXBuilder.append("sex=1_B").append("\t");
            headerMixedXBuilder.append("sex=1_AA").append("\t");
            headerMixedXBuilder.append("sex=1_AB").append("\t");
            headerMixedXBuilder.append("sex=1_BB").append("\t");
            headerMixedXBuilder.append("sex=1_NULL").append("\t");
            headerMixedXBuilder.append("sex=1_total").append("\t");
            headerMixedXBuilder.append("sex=1_maf").append("\t");
            headerMixedXBuilder.append("sex=1_info").append("\t");
            headerMixedXBuilder.append("sex=1_impute_info").append("\t");
        } else {
            headerMixedXBuilder.append("sex=2_A").append("\t");
            headerMixedXBuilder.append("sex=2_B").append("\t");
            headerMixedXBuilder.append("sex=2_AA").append("\t");
            headerMixedXBuilder.append("sex=2_AB").append("\t");
            headerMixedXBuilder.append("sex=2_BB").append("\t");
            headerMixedXBuilder.append("sex=2_NULL").append("\t");
            headerMixedXBuilder.append("sex=2_total").append("\t");
            headerMixedXBuilder.append("sex=2_maf").append("\t");
            headerMixedXBuilder.append("sex=2_info").append("\t");
            headerMixedXBuilder.append("sex=2_impute_info").append("\t");	
        }
        headerMixedXBuilder.append("frequentist_add_null_ll").append("\t");
        headerMixedXBuilder.append("frequentist_add_alternative_ll").append("\t");
        if(GuidanceImpl.getSex1().equals(sex)) {
        	headerMixedXBuilder.append("frequentist_add_beta_1:genotype/sex=1").append("\t");
        	headerMixedXBuilder.append("frequentist_add_se_1:genotype/sex=1").append("\t");
        } else {
            headerMixedXBuilder.append("frequentist_add_beta_2:genotype/sex=2").append("\t");
            headerMixedXBuilder.append("frequentist_add_se_2:genotype/sex=2").append("\t");	
        }
        headerMixedXBuilder.append("frequentist_add_degrees_of_freedom").append("\t");
        headerMixedXBuilder.append("frequentist_add_pvalue").append("\t");
        headerMixedXBuilder.append("comment");

        return headerMixedXBuilder.toString();
    }

    /**
     * Construct the header (1-22)
     * 
     * @return
     */
    public static String constructHeader() {
        StringBuilder headerMixedBuilder = new StringBuilder();
        headerMixedBuilder.append("chr").append("\t");
        headerMixedBuilder.append("position").append("\t");
        headerMixedBuilder.append("rs_id_all").append("\t");
        headerMixedBuilder.append("info_all").append("\t");
        headerMixedBuilder.append("certainty_all").append("\t");
        headerMixedBuilder.append("alleleA").append("\t");
        headerMixedBuilder.append("alleleB").append("\t");
        headerMixedBuilder.append("index").append("\t");
        headerMixedBuilder.append("average_maximum_posterior_call").append("\t");
        headerMixedBuilder.append("info").append("\t");
        headerMixedBuilder.append("cohort_1_AA").append("\t");
        headerMixedBuilder.append("cohort_1_AB").append("\t");
        headerMixedBuilder.append("cohort_1_BB").append("\t");
        headerMixedBuilder.append("cohort_1_NULL").append("\t");
        headerMixedBuilder.append("all_AA").append("\t");
        headerMixedBuilder.append("all_AB").append("\t");
        headerMixedBuilder.append("all_BB").append("\t");
        headerMixedBuilder.append("all_NULL").append("\t");
        headerMixedBuilder.append("all_total").append("\t");
        headerMixedBuilder.append("cases_AA").append("\t");
        headerMixedBuilder.append("cases_AB").append("\t");
        headerMixedBuilder.append("cases_BB").append("\t");
        headerMixedBuilder.append("cases_NULL").append("\t");
        headerMixedBuilder.append("cases_total").append("\t");
        headerMixedBuilder.append("controls_AA").append("\t");
        headerMixedBuilder.append("controls_AB").append("\t");
        headerMixedBuilder.append("controls_BB").append("\t");
        headerMixedBuilder.append("controls_NULL").append("\t");
        headerMixedBuilder.append("controls_total").append("\t");
        headerMixedBuilder.append("all_maf").append("\t");
        headerMixedBuilder.append("cases_maf").append("\t");
        headerMixedBuilder.append("controls_maf").append("\t");
        headerMixedBuilder.append("missing_data_proportion").append("\t");
        headerMixedBuilder.append("cohort_1_hwe").append("\t");
        headerMixedBuilder.append("cases_hwe").append("\t");
        headerMixedBuilder.append("controls_hwe").append("\t");
        headerMixedBuilder.append("het_OR").append("\t");
        headerMixedBuilder.append("het_OR_lower").append("\t");
        headerMixedBuilder.append("het_OR_upper").append("\t");
        headerMixedBuilder.append("hom_OR").append("\t");
        headerMixedBuilder.append("hom_OR_lower").append("\t");
        headerMixedBuilder.append("hom_OR_upper").append("\t");
        headerMixedBuilder.append("all_OR").append("\t");
        headerMixedBuilder.append("all_OR_lower").append("\t");
        headerMixedBuilder.append("all_OR_upper").append("\t");
        headerMixedBuilder.append("frequentist_add_pvalue").append("\t");
        headerMixedBuilder.append("frequentist_add_info").append("\t");
        headerMixedBuilder.append("frequentist_add_beta_1").append("\t");
        headerMixedBuilder.append("frequentist_add_se_1").append("\t");
        headerMixedBuilder.append("frequentist_dom_pvalue").append("\t");
        headerMixedBuilder.append("frequentist_dom_info").append("\t");
        headerMixedBuilder.append("frequentist_dom_beta_1").append("\t");
        headerMixedBuilder.append("frequentist_dom_se_1").append("\t");
        headerMixedBuilder.append("frequentist_rec_pvalue").append("\t");
        headerMixedBuilder.append("frequentist_rec_info").append("\t");
        headerMixedBuilder.append("frequentist_rec_beta_1").append("\t");
        headerMixedBuilder.append("frequentist_rec_se_1").append("\t");
        headerMixedBuilder.append("frequentist_gen_pvalue").append("\t");
        headerMixedBuilder.append("frequentist_gen_info").append("\t");
        headerMixedBuilder.append("frequentist_gen_beta_1").append("\t");
        headerMixedBuilder.append("frequentist_gen_se_1").append("\t");
        headerMixedBuilder.append("frequentist_gen_beta_2").append("\t");
        headerMixedBuilder.append("frequentist_gen_se_2").append("\t");
        headerMixedBuilder.append("frequentist_het_pvalue").append("\t");
        headerMixedBuilder.append("frequentist_het_info").append("\t");
        headerMixedBuilder.append("frequentist_het_beta_1").append("\t");
        headerMixedBuilder.append("frequentist_het_se_1").append("\t");
        headerMixedBuilder.append("comment");
        return headerMixedBuilder.toString();
    }

    /**
     * Construct the condensed header
     * 
     * @return
     */
    public static String constructCondensedHeader() {
        StringBuilder condensedHeaderBuilder = new StringBuilder();
        condensedHeaderBuilder.append("chr").append("\t");
        condensedHeaderBuilder.append("position").append("\t");
        condensedHeaderBuilder.append("alleleA").append("\t");
        condensedHeaderBuilder.append("alleleB").append("\t");
        condensedHeaderBuilder.append("frequentist_add_pvalue").append("\t");
        condensedHeaderBuilder.append("info_all");
        return condensedHeaderBuilder.toString();
    }

    /**
     * Method to create a HashMap from the header of particular files -> from header to index
     * 
     * @param line
     * @param separator
     * @return
     */
    public static HashMap<String, Integer> createHashWithHeader(String line, String separator) {
        HashMap<String, Integer> myHashLine = new HashMap<>();
        System.out.println("Starting hash creation");
        String[] splitted = line.split(separator);
        for (int i = 0; i < splitted.length; i++) {
        	System.out.println("Index: " + Integer.toString(i) + " has value " + splitted[i]);
            myHashLine.put(splitted[i], i);
        }
        return myHashLine;
    }

    /**
     * Method to create a HashMap from the header of particular files -> from index to header
     * 
     * @param line
     * @param separator
     * @return
     */
    public static HashMap<Integer, String> createHashWithHeaderReversed(String line, String separator) {
        HashMap<Integer, String> myHashLine = new HashMap<>();

        String[] splitted = line.split(separator);
        for (int i = 0; i < splitted.length; i++) {
            myHashLine.put(i, splitted[i]);
        }
        return myHashLine;
    }

}
