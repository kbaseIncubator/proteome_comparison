package genomeproteomecomparison;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map;
import us.kbase.auth.AuthToken;
import us.kbase.common.service.JsonServerMethod;
import us.kbase.common.service.JsonServerServlet;
import us.kbase.common.service.JsonServerSyslog;
import us.kbase.common.service.RpcContext;

//BEGIN_HEADER
//END_HEADER

/**
 * <p>Original spec-file module name: GenomeProteomeComparison</p>
 * <pre>
 * </pre>
 */
public class GenomeProteomeComparisonServer extends JsonServerServlet {
    private static final long serialVersionUID = 1L;
    private static final String version = "0.0.8";
    private static final String gitUrl = "https://github.com/mrcreosote/genomeproteomecomparison.git";
    private static final String gitCommitHash = "ebf1ab51829c5ed2abd3c793b4f3637e6feb503d";

    //BEGIN_CLASS_HEADER
    //END_CLASS_HEADER

    public GenomeProteomeComparisonServer() throws Exception {
        super("GenomeProteomeComparison");
        //BEGIN_CONSTRUCTOR
        //END_CONSTRUCTOR
    }

    /**
     * <p>Original spec-file function name: blast_proteomes</p>
     * <pre>
     * </pre>
     * @param   input   instance of type {@link genomeproteomecomparison.BlastProteomesParams BlastProteomesParams} (original type "blast_proteomes_params")
     * @return   parameter "output_ref" of original type "ws_protcmp_id" (A workspace ID that references a Genome data object. @id ws ProteomeComparison)
     */
    @JsonServerMethod(rpc = "GenomeProteomeComparison.blast_proteomes", async=true)
    public String blastProteomes(BlastProteomesParams input, AuthToken authPart, RpcContext jsonRpcContext) throws Exception {
        String returnVal = null;
        //BEGIN blast_proteomes
        returnVal = new BlastProteomes(authPart, config).run(input);
        //END blast_proteomes
        return returnVal;
    }
    @JsonServerMethod(rpc = "GenomeProteomeComparison.status")
    public Map<String, Object> status() {
        Map<String, Object> returnVal = null;
        //BEGIN_STATUS
        returnVal = new LinkedHashMap<String, Object>();
        returnVal.put("state", "OK");
        returnVal.put("message", "");
        returnVal.put("version", version);
        returnVal.put("git_url", gitUrl);
        returnVal.put("git_commit_hash", gitCommitHash);
        //END_STATUS
        return returnVal;
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 1) {
            new GenomeProteomeComparisonServer().startupServer(Integer.parseInt(args[0]));
        } else if (args.length == 3) {
            JsonServerSyslog.setStaticUseSyslog(false);
            JsonServerSyslog.setStaticMlogFile(args[1] + ".log");
            new GenomeProteomeComparisonServer().processRpcCall(new File(args[0]), new File(args[1]), args[2]);
        } else {
            System.out.println("Usage: <program> <server_port>");
            System.out.println("   or: <program> <context_json_file> <output_json_file> <token>");
            return;
        }
    }
}
