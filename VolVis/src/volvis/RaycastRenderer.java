/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;



/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    static int MAXINT = 1000000;
    Mode mode = Mode.slicer;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    
    private VoxelGradient getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] > volume.getDimX()-1 || coord[1] < 0 || coord[1] > volume.getDimY()-1
                || coord[2] < 0 || coord[2] > volume.getDimZ()-1) {
            return new VoxelGradient();
        }
        return gradients.getGradient((int)coord[0], (int)coord[1], (int)coord[2]);
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] <= 0 || coord[0] >= volume.getDimX()-1 || coord[1] <= 0 || coord[1] >= volume.getDimY()-1
                || coord[2] <= 0 || coord[2] >= volume.getDimZ()-1) {
            return 0;
        }
        double alpha = coord[0] - Math.floor(coord[0]);
        double beta = coord[1] - Math.floor(coord[1]);
        double gamma = coord[2] - Math.floor(coord[2]);
        
        int x0 = volume.getVoxel((int) Math.floor(coord[0]), (int) Math.floor(coord[1]), (int) Math.floor(coord[2]));
        int x1 = volume.getVoxel((int) Math.ceil(coord[0]), (int) Math.floor(coord[1]), (int) Math.floor(coord[2]));
        int x2 = volume.getVoxel((int) Math.floor(coord[0]), (int) Math.ceil(coord[1]), (int) Math.floor(coord[2]));
        int x3 = volume.getVoxel((int) Math.ceil(coord[0]), (int) Math.ceil(coord[1]), (int) Math.floor(coord[2]));
        int x4 = volume.getVoxel((int) Math.floor(coord[0]), (int) Math.floor(coord[1]), (int) Math.ceil(coord[2]));
        int x5 = volume.getVoxel((int) Math.ceil(coord[0]), (int) Math.floor(coord[1]), (int) Math.ceil(coord[2]));
        int x6 = volume.getVoxel((int) Math.floor(coord[0]), (int) Math.ceil(coord[1]), (int) Math.ceil(coord[2]));
        int x7 = volume.getVoxel((int) Math.ceil(coord[0]), (int) Math.ceil(coord[1]), (int) Math.ceil(coord[2]));
        
        double result = (1.0-alpha)*(1.0-beta)*(1.0-gamma)*(double)x0 + alpha*(1.0-beta)*(1.0-gamma)*(double)x1 
                + (1.0-alpha)*beta*(1.0-gamma)*(double)x2 + alpha*beta*(1.0-gamma)*(double)x3
                + (1.0-alpha)*(1.0-beta)*gamma*(double)x4 + alpha*(1.0-beta)*gamma*(double)x5
                +(1.0-alpha)*beta*gamma*(double)x6 + alpha*beta*gamma*(double)x6;
        return (short) result; //NAVRAGEN!! Casting?
    }
    

    void calculate(double[] viewMatrix, Mode mode) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        if (mode == Mode.slicer){
            for (int j = 0; j < image.getHeight(); j++) {
                for (int i = 0; i < image.getWidth(); i++) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2];           
                int val = getVoxel(pixelCoord);
                voxelColor = tFunc.getColor(val);
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                }
            }
        }
        else{
            for (int j = 0; j < image.getHeight(); j++) {
                for (int i = 0; i < image.getWidth(); i++) {
                    int maxDim = Math.max(Math.max(volume.getDimX(), volume.getDimY()), volume.getDimZ());
                    short [] values = new short[maxDim];
                    VoxelGradient [] grads = new VoxelGradient[maxDim];
                    for (int k = 0; k < maxDim; k = k+1){
                        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k - imageCenter)
                                + volumeCenter[0];
                        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k - imageCenter)
                                + volumeCenter[1];
                        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k - imageCenter)
                                + volumeCenter[2];
                        grads[k] = getGradient(pixelCoord);
                        values[k] = getVoxel(pixelCoord);
                    }
                if (mode == Mode.mip){
                    int val = mip(values);
                    // Map the intensity to a grey value by linear scaling
                    voxelColor.r = val/max;
                    voxelColor.g = voxelColor.r;
                    voxelColor.b = voxelColor.r;
                    voxelColor.a = val > 0 ? 1.0 : 0.0;
                }
                else if (mode == Mode.compositing)
                    voxelColor = compositing(values);
                else if (mode == Mode.transfer2d)
                    voxelColor = opacityCalc(values, grads, viewVec);
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                }
            }
        }
    }
    
    private int mip(short [] values){
        short max = 0;
        for (int i = 0; i < values.length; i++){
            if(values[i]>max)
                max = values[i];
        }
        return max;
    }
    
    private TFColor compositing(short [] values){
        TFColor [] colors = new TFColor[values.length];
        for(int i =0; i < values.length; i++)
            colors[i] = tFunc.getColor(values[i]);
        
        double [] startvalues = new double[]{0.0, 0.0,0.0};
        double [] compvalues = comp(colors, colors.length-1, startvalues); 
        
        return new TFColor(compvalues[0], compvalues[1], compvalues[2], 1.0);
    }
    
    private double [] comp(TFColor[] colors, int i, double [] oldvalue){
        if(i<=0)
            return oldvalue;
        else{
            double [] value = new double [3];
            value[0] = colors[i].r*colors[i].a + (1-colors[i].a)*oldvalue[0];
            value[1] = colors[i].g*colors[i].a + (1-colors[i].a)*oldvalue[1];
            value[2] = colors[i].b*colors[i].a + (1-colors[i].a)*oldvalue[2];
            return comp(colors, (i-1), value);
        }
    }
    
    private TFColor opacityCalc(short [] values, VoxelGradient [] gradients, double[] viewVec){
        TFColor c = tfEditor2D.triangleWidget.color;
        double r = tfEditor2D.triangleWidget.radius;
        short fv = tfEditor2D.triangleWidget.baseIntensity;
        TFColor [] colors = new TFColor[values.length];
        
        
        for (int i =0; i<values.length; i++){
            float gradmag = gradients[i].mag;
            c = colorCalc(c,  gradients[i], viewVec);
            if (gradmag == 0 && values[i] == fv)
                colors[i] = new TFColor(c.r,c.g,c.b,c.a*1);
            else if (gradients[i].mag > 0 && (values[i]-(r*gradmag)<= fv && values[i]+(r*gradmag) >= fv)){
                double alpha = 1 - (1/r)*((fv-values[i])/gradmag);
                colors[i] = new TFColor(c.r,c.g,c.b,c.a*alpha);
            }
            else
                colors[i] = new TFColor(c.r,c.g,c.b,c.a*0);
        }
        double [] startvalues = new double[]{1.0, 1.0,1.0};
        double [] compvalues = comp(colors, colors.length-1, startvalues); 
        
        return new TFColor(compvalues[0], compvalues[1], compvalues[2], 1.0);
    }
    
    private TFColor colorCalc(TFColor c, VoxelGradient grad, double[] viewVec){
        double magV = Math.sqrt(c.r*c.r+c.g*c.g+c.b*c.b);
        double [] H = new double[]{(2*viewVec[0])/magV, (2*viewVec[1])/magV, (2*viewVec[2])/magV};
        float [] N = new float[]{grad.x/grad.mag, grad.y/grad.mag, grad.z/grad.mag};
        double Ia = 0.1;
        double kdiff = 0.7;
        double kspec = 0.2;
        double alpha = 10;
        double r = 0.0;//Ia + c.r*kdiff*();
        double g = 0.0;
        double b = 0.0;
        return new TFColor(r,g,b,c.a);        
    }
    
    
    
//    private TFColor Levoy (short [] values){
//        TFColor [] colors = new TFColor[values.length];
//        for(int i =0; i < values.length; i++)
//            colors[i] = tFunc.getColor(values[i]);
//        
//        return new TFColor();
//    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }
    
    private enum Mode {
    slicer, mip, compositing, transfer2d;
    }
    
    public void setMode(String mode){
        Mode m = Mode.valueOf(mode);
        this.mode = m;
    }

    @Override
    public void visualize(GL2 gl) {
        //TODO: Change this:
        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        calculate(viewMatrix, mode);          
           
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
