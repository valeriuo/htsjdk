package htsjdk.variant.vcf;

import htsjdk.HtsjdkTest;
import htsjdk.samtools.util.IOUtil;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Paths;

/**
 * Created by farjoun on 10/12/17.
 */
public class VCFPathReaderTest extends HtsjdkTest {

    @DataProvider(name = "pathsData")
    Object[][] pathsData() {
        return new Object[][]{
                // various ways to refer to a local file
                {"src/test/resources/htsjdk/variant/VCF4HeaderTest.vcf", false, true},

                {Paths.get("").toAbsolutePath().toString() + "/src/test/resources/htsjdk/variant/VCF4HeaderTest.vcf", false, true},
                {"file://" + Paths.get("").toAbsolutePath().toString() + "/src/test/resources/htsjdk/variant/VCF4HeaderTest.vcf", false, true},

                //testing GCS files:

                // this is almost a vcf, but not quite it's missing the #CHROM line and it has no content...
                {"gs://htsjdk-testdata/htsjdk/variant/Homo_sapiens_assembly38.tile_db_header.vcf", false, false},

                // test that have indexes
                {"gs://htsjdk-testdata/htsjdk/variant/test.vcf.bgz", true, true},
                {"gs://htsjdk-testdata/htsjdk/variant/serialization_test.bcf", true, true},
                {"gs://htsjdk-testdata/htsjdk/variant/test1.vcf", true, true},

                // test that lack indexes (should succeed)
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.vcf.gz", false, true},
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.vcf", false, true},
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.bcf", false, true},
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.vcf.bgz", false, true},

                // test that lack indexes should fail)
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.vcf.gz", true, false},
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.vcf", true, false},
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.bcf", true, false},
                {"gs://htsjdk-testdata/htsjdk/variant/VcfThatLacksAnIndex.vcf.bgz", true, false},

                // testing a non-existent scheme:
                {"bogus://" + Paths.get("").toAbsolutePath().toString() + "/src/test/resources/htsjdk/variant/VCF4HeaderTest.vcf", false, false},
        };
    }

    @Test(dataProvider = "pathsData", timeOut = 60_000)
    public void testCanOpenVCFPathReader(final String uri, final boolean requiresIndex, final boolean shouldSucceed) throws Exception {

        // read an existing VCF
        try (final VCFPathReader reader = new VCFPathReader(IOUtil.getPath(uri), requiresIndex)) {
            final VCFHeader header = reader.getFileHeader();
        } catch (Exception e) {
            if (shouldSucceed) throw e;
        }
    }
}