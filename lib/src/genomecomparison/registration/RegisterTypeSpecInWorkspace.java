package genomecomparison.registration;

import genomecomparison.Utils;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.net.URL;

import us.kbase.auth.AuthService;
import us.kbase.workspace.RegisterTypespecParams;
import us.kbase.workspace.WorkspaceClient;

public class RegisterTypeSpecInWorkspace {
    private static final String wsUrl = "https://ci.kbase.us/services/ws";
    private static final String moduleName = "GenomeComparison";
    private static final String owner = "rsutormin";
    
    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            System.err.println("Usage: <program> <owner-password>");
            System.exit(1);
        }
        WorkspaceClient wc = new WorkspaceClient(new URL(wsUrl),
                AuthService.login(owner, args[0]).getToken());
        wc.setIsInsecureHttpConnectionAllowed(true);
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        Utils.copy(new FileInputStream(new File(moduleName + ".spec")), baos);
        String spec = new String(baos.toByteArray());
        System.out.println(wc.registerTypespec(new RegisterTypespecParams().withSpec(spec)
                .withDryrun(0L)));
        System.out.println(wc.releaseModule(moduleName));
    }
}
