package genomeproteomecomparison.test;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import junit.framework.Assert;

import org.ini4j.Ini;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import genomecomparison.ProteomeComparison;
import genomeproteomecomparison.BlastProteomesParams;
import genomeproteomecomparison.GenomeProteomeComparisonServer;
import us.kbase.auth.AuthService;
import us.kbase.auth.AuthToken;
import us.kbase.common.service.JsonServerSyslog;
import us.kbase.common.service.RpcContext;
import us.kbase.common.service.UObject;
import us.kbase.workspace.CreateWorkspaceParams;
import us.kbase.workspace.GetObjects2Params;
import us.kbase.workspace.ObjectSaveData;
import us.kbase.workspace.ObjectSpecification;
import us.kbase.workspace.ProvenanceAction;
import us.kbase.workspace.SaveObjectsParams;
import us.kbase.workspace.WorkspaceClient;
import us.kbase.workspace.WorkspaceIdentity;

public class GenomeComparisonServerTest {
    private static AuthToken token = null;
    private static Map<String, String> config = null;
    private static WorkspaceClient wsClient = null;
    private static String wsName = null;
    private static GenomeProteomeComparisonServer impl = null;
    private static String genome1objName = "Ecoli_042.genome";
    private static String genome2objName = "Ecoli_K12.genome";
    
    @BeforeClass
    public static void init() throws Exception {
        token = AuthService.validateToken(System.getenv("KB_AUTH_TOKEN"));
        String configFilePath = System.getenv("KB_DEPLOYMENT_CONFIG");
        File deploy = new File(configFilePath);
        Ini ini = new Ini(deploy);
        config = ini.get("GenomeProteomeComparison");
        wsClient = new WorkspaceClient(new URL(config.get("workspace-url")), token);
        wsClient.setIsInsecureHttpConnectionAllowed(true);
        // These lines are necessary because we don't want to start linux syslog bridge service
        JsonServerSyslog.setStaticUseSyslog(false);
        JsonServerSyslog.setStaticMlogFile(new File(config.get("scratch"), "test.log").getAbsolutePath());
        impl = new GenomeProteomeComparisonServer();
        // Upload genomes
        String[] genomeObjNames = {genome1objName, genome2objName};
        for (String genomeObjName : genomeObjNames) {
            String contigsetObjName = genomeObjName + ".contigset";
            Map<String, Object> contigsetData = new LinkedHashMap<String, Object>();
            contigsetData.put("contigs", new ArrayList<Object>());
            contigsetData.put("id", contigsetObjName);
            contigsetData.put("md5", "md5");
            contigsetData.put("name", contigsetObjName);
            contigsetData.put("source", "User uploaded data");
            contigsetData.put("source_id", "noid");
            contigsetData.put("type", "Organism");
            wsClient.saveObjects(new SaveObjectsParams().withWorkspace(getWsName()).withObjects(Arrays.asList(
                    new ObjectSaveData().withName(contigsetObjName).withType("KBaseGenomes.ContigSet")
                    .withData(new UObject(contigsetData)))));
            File inputDir = new File("test/data");
            File inputFile = new File(inputDir, genomeObjName + ".json.gz");
            InputStream is = new GZIPInputStream(new FileInputStream(inputFile));
            @SuppressWarnings("unchecked")
            Map<String, Object> genomeData = UObject.getMapper().readValue(is, Map.class);
            is.close();
            genomeData.put("contigset_ref", getWsName() + "/" + contigsetObjName);
            wsClient.saveObjects(new SaveObjectsParams().withWorkspace(getWsName()).withObjects(Arrays.asList(
                    new ObjectSaveData().withName(genomeObjName).withType("KBaseGenomes.Genome")
                    .withData(new UObject(genomeData)))));
        }
    }
    
    private static String getWsName() throws Exception {
        if (wsName == null) {
            long suffix = System.currentTimeMillis();
            wsName = "test_GenomeProteomeComparison_" + suffix;
            wsClient.createWorkspace(new CreateWorkspaceParams().withWorkspace(wsName));
        }
        return wsName;
    }
    
    private static RpcContext getContext() {
        return new RpcContext().withProvenance(Arrays.asList(new ProvenanceAction()
            .withService("GenomeProteomeComparison").withMethod("please_never_use_it_in_production")
            .withMethodParams(new ArrayList<UObject>())));
    }
    
    @AfterClass
    public static void cleanup() {
        if (wsName != null) {
            try {
                wsClient.deleteWorkspace(new WorkspaceIdentity().withWorkspace(wsName));
                System.out.println("Test workspace was deleted");
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }
    
    @Test
    public void testBlastProteomesOk() throws Exception {
        String outputObjName = "protcmp.1";
        impl.blastProteomes(
                new BlastProteomesParams().withGenome1id(genome1objName).withGenome1ws(getWsName())
                .withGenome2id(genome2objName).withGenome2ws(getWsName()).withOutputId(outputObjName)
                .withOutputWs(getWsName()), token, getContext());
        ProteomeComparison protcmp = wsClient.getObjects2(new GetObjects2Params().withObjects(Arrays.asList(
                new ObjectSpecification().withWorkspace(getWsName()).withName(outputObjName))))
                .getData().get(0).getData().asClassInstance(ProteomeComparison.class);
        int size = 0;
        for (List<?> row : protcmp.getData1())
            size += row.size();
        Assert.assertTrue(size > 4000);
    }
}