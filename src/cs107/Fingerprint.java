package cs107;

import javax.crypto.spec.PSource;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cs107.Main.arrayEqual;
import static cs107.Main.printArray;

/**
 * Provides tools to compare fingerprint.
 */
public class Fingerprint {

    /**
     * The number of pixels to consider in each direction when doing the linear
     * regression to compute the orientation.
     */
    public static final int ORIENTATION_DISTANCE = 16;

    /**
     * The maximum distance between two minutiae to be considered matching.
     */
    public static final int DISTANCE_THRESHOLD = 5;

    /**
     * The number of matching minutiae needed for two fingerprints to be considered
     * identical.
     */
    public static final int FOUND_THRESHOLD = 20;

    /**
     * The distance between two angle to be considered identical.
     */
    public static final int ORIENTATION_THRESHOLD = 20;

    /**
     * The offset in each direction for the rotation to test when doing the
     * matching.
     */
    public static final int MATCH_ANGLE_OFFSET = 2;

    /**
     * Returns an array containing the value of the 8 neighbours of the pixel at
     * coordinates <code>(row, col)</code>.
     * <p>
     * The pixels are returned such that their indices corresponds to the following
     * diagram:<br>
     * ------------- <br>
     * | 7 | 0 | 1 | <br>
     * ------------- <br>
     * | 6 | _ | 2 | <br>
     * ------------- <br>
     * | 5 | 4 | 3 | <br>
     * ------------- <br>
     * <p>
     * If a neighbours is out of bounds of the image, it is considered white.
     * <p>
     * If the <code>row</code> or the <code>col</code> is out of bounds of the
     * image, the returned value should be <code>null</code>.
     *
     * @param image array containing each pixel's boolean value.
     * @param row   the row of the pixel of interest, must be between
     *              <code>0</code>(included) and
     *              <code>image.length</code>(excluded).
     * @param col   the column of the pixel of interest, must be between
     *              <code>0</code>(included) and
     *              <code>image[row].length</code>(excluded).
     * @return An array containing each neighbours' value.
     */

    // Function findBoolean: tests whether each pixel is in-bounds
    public static boolean findBoolean(boolean[][] image, int row, int col) {

        return (row >= 0) && (row <= image.length - 1) && (col >= 0)
                && (col <= image[row].length - 1) && image[row][col];
    }

    public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
        assert (image != null); // special case that is not expected (the image is supposed to have been checked
        // earlier)
        //TODO implement
        boolean[] neighbours = new boolean[8];

        neighbours[0] = findBoolean(image,row-1,col);
        neighbours[1] = findBoolean(image, row-1, col+1);
        neighbours[2] = findBoolean(image, row,col+1);
        neighbours[3] = findBoolean(image, row+1,col+1);
        neighbours[4] = findBoolean(image, row+1,col);
        neighbours[5] = findBoolean(image,row+1,col-1);
        neighbours[6] = findBoolean(image,row,col-1);
        neighbours[7] = findBoolean(image,row-1,col-1);

        return neighbours;
    }

    /**
     * Computes the number of black (<code>true</code>) pixels among the neighbours
     * of a pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @return the number of black neighbours.
     */
    public static int blackNeighbours(boolean[] neighbours) {
        //TODO implement

        int numberOfBlackNeighbours = 0;
        for (boolean blackPixel : neighbours) if (blackPixel) numberOfBlackNeighbours++;
        return numberOfBlackNeighbours;
    }

    /**
     * Computes the number of white to black transitions among the neighbours of
     * pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @return the number of white to black transitions.
     */
    public static int transitions(boolean[] neighbours) {
        //TODO implement
        int transitions = 0;
        for (int i = 0 ; i < 7 ; i++) {
            if (!neighbours[i] && neighbours[i+1]) {
                transitions++;
            }
        }
        if (!neighbours[7] && neighbours[0]) {
            transitions++;
        }
        return transitions;
    }

    /**
     * Returns <code>true</code> if the images are identical and false otherwise.
     *
     * @param image1 array containing each pixel's boolean value.
     * @param image2 array containing each pixel's boolean value.
     * @return <code>True</code> if they are identical, <code>false</code>
     *         otherwise.
     */
    public static boolean identical(boolean[][] image1, boolean[][] image2) {
        //TODO implement
        return arrayEqual(image1, image2);
    }

    //added something
    /**
     * Internal method used by {@link #thin(boolean[][])}.
     *
     * @param image array containing each pixel's boolean value.
     * @param step  the step to apply, Step 0 or Step 1.
     * @return A new array containing each pixel's value after the step.
     */
    public static boolean[][] thinningStep(boolean[][] image, int step) {
        //TODO implement

        // Working with copy of image
        boolean[][] copy = copy(image);

        if (step == 0) {
            for (int i=0;i<image.length;++i) {
                for (int j=0;j<image[i].length;++j) {
                    boolean[] neighbours = getNeighbours(image, i, j);
                    if (image[i][j] && (blackNeighbours(neighbours)<=6)
                            && (blackNeighbours(neighbours)>=2) && (transitions(neighbours)==1)
                            && !(neighbours[0] && neighbours[2] && neighbours[4])
                            && !(neighbours[2] && neighbours[4] && neighbours[6])) {
                        copy[i][j] = false;
                    }
                }
            }
        } else if (step == 1) {
            for (int i=0;i<image.length;++i) {
                for (int j=0;j<image[i].length;++j) {
                    boolean[] neighbours = getNeighbours(image, i, j);
                    if (image[i][j] && (blackNeighbours(neighbours)<=6)
                            && (blackNeighbours(neighbours)>=2) && (transitions(neighbours)==1)
                            && !(neighbours[0] && neighbours[2] && neighbours[6])
                            && !(neighbours[0] && neighbours[4] && neighbours[6])) {
                        copy[i][j] = false;
                    }
                }
            }
        }
        return copy;
    }

    // Function copy: makes a copy of an image
    public static boolean[][] copy(boolean[][]image) {
        boolean[][] copy = new boolean[image.length][image[0].length];
        for (int i = 0; i<image.length; ++i) {
            for (int j =0; j<image[0].length; ++j) {
                copy[i][j] = image[i][j];
            }
        }
        return copy;
    }

    /**
     * Compute the skeleton of a boolean image.
     *
     * @param image array containing each pixel's boolean value.
     * @return array containing the boolean value of each pixel of the image after
     *         applying the thinning algorithm.
     */
    public static boolean[][] thin(boolean[][] image) {
        //TODO implement
        boolean[][]imageDepart;

        do {
            imageDepart = copy(image);
            image = thinningStep(image,0);
            image = thinningStep(image,1);
        } while (!identical(imageDepart,image));

        return image;
    }

    /**
     * Computes all pixels that are connected to the pixel at coordinate
     * <code>(row, col)</code> and within the given distance of the pixel.
     *
     * @param image    array containing each pixel's boolean value.
     * @param row      the first coordinate of the pixel of interest.
     * @param col      the second coordinate of the pixel of interest.
     * @param distance the maximum distance at which a pixel is considered.
     * @return An array where <code>true</code> means that the pixel is within
     *         <code>distance</code> and connected to the pixel at
     *         <code>(row, col)</code>.
     */


    public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {
        //TODO implement

        boolean[][]pixRes = new boolean[image.length][image[0].length];
        pixRes[row][col] = true;
        boolean pixChanged;

        do {
            pixChanged = false;
            for (int i = 0 ; i < image.length ; i++) {
                for (int j = 0 ; j < image[0].length ; j++) {
                    if ((image[i][j]) && (blackNeighbours(getNeighbours(pixRes,i,j)) > 0)
                            && !(i < row - distance || i > row + distance
                              || j < col - distance || j > col + distance)) {

                        if (!pixRes[i][j]) {
                            pixRes[i][j] = true;
                            pixChanged = true;
                        }
                    }
                }
            }
        } while (pixChanged);
        return pixRes;
    }

    /**
     * Computes the slope of a minutia using linear regression.
     *
     * @param connectedPixels the result of
     *                        {@link #connectedPixels(boolean[][], int, int, int)}.
     * @param row             the row of the minutia.
     * @param col             the col of the minutia.
     * @return the slope.
     */
    public static double computeSlope(boolean[][] connectedPixels, int row, int col) {
        //TODO implement
        double sommeX = 0;
        double sommeY = 0;
        double sommeXY = 0;

        for (int i = 0; i < connectedPixels.length; ++i) {
            for (int j = 0; j < connectedPixels[0].length; ++j) {
                if (connectedPixels[i][j]) {
                    sommeX += Math.pow(j - col, 2);
                    sommeY += Math.pow(row - i, 2);
                    sommeXY += (j - col) * (row - i);
                }
            }
        }

        return (sommeX == 0) ? Double.POSITIVE_INFINITY : (sommeX < sommeY) ? (sommeY / sommeXY) : (sommeXY / sommeX);
    }


    /**
     * Computes the orientation of a minutia in radians.
     *
     * @param connectedPixels the result of
     *                        {@link #connectedPixels(boolean[][], int, int, int)}.
     * @param row             the row of the minutia.
     * @param col             the col of the minutia.
     * @param slope           the slope as returned by
     *                        {@link #computeSlope(boolean[][], int, int)}.
     * @return the orientation of the minutia in radians.
     */
    public static double computeAngle(boolean[][] connectedPixels, int row, int col, double slope) {
        //TODO implement

        double angle = Math.atan(slope);
        int dessus = 0;
        int dessous = 0;

        for (int i = 0; i < connectedPixels.length; ++i) {
            for (int j = 0; j < connectedPixels[0].length; ++j) {
                if (connectedPixels[i][j]) {
                    if (row - i >= -(j - col) / slope) {
                        ++dessus;
                    } else {
                        ++dessous;
                    }
                }
            }
        }
        if (((dessous > dessus) && (angle > 0)) || ((angle < 0) && (dessous < dessus))) {angle += Math.PI;}
        else if (slope == Double.POSITIVE_INFINITY) {angle = (dessous < dessus) ? Math.PI/2 : -Math.PI/2;}

        return angle;
    }

    /**
     * Computes the orientation of the minutia that the coordinate <code>(row,
     * col)</code>.
     *
     * @param image    array containing each pixel's boolean value.
     * @param row      the first coordinate of the pixel of interest.
     * @param col      the second coordinate of the pixel of interest.
     * @param distance the distance to be considered in each direction to compute
     *                 the orientation.
     * @return The orientation in degrees.
     */
    public static int computeOrientation(boolean[][] image, int row, int col, int distance) {
        //TODO implement
        boolean[][] pixelsProches = connectedPixels(image,row,col,ORIENTATION_DISTANCE);
        double pente = computeSlope(pixelsProches, row, col);
        double angle = computeAngle(pixelsProches, row, col, pente);

        angle = Math.round(Math.toDegrees(angle));

        if (angle < 0) {angle = angle + 360;}

        return (int)angle;
    }

    /**
     * Extracts the minutiae from a thinned image.
     *
     * @param image array containing each pixel's boolean value.
     * @return The list of all minutiae. A minutia is represented by an array where
     *         the first element is the row, the second is column, and the third is
     *         the angle in degrees.
     * @see #thin(boolean[][])
     */
    public static List<int[]> extract(boolean[][] image) {
        //TODO implement
        List<int[]> minutiae = new ArrayList<>();

        for (int i = 1; i < image.length-1; ++i) {
            for (int j = 1; j < image[0].length-1; ++j) {
                if (image[i][j]
                        && (transitions(getNeighbours(image, i, j)) == 3
                        || (transitions(getNeighbours(image, i, j)) == 1))) {
                    minutiae.add(new int[] {i,j,computeOrientation(image,i,j,ORIENTATION_DISTANCE)});
                }
            }
        }

        return minutiae;
    }

    /**
     * Applies the specified rotation to the minutia.
     *
     * @param minutia   the original minutia.
     * @param centerRow the row of the center of rotation.
     * @param centerCol the col of the center of rotation.
     * @param rotation  the rotation in degrees.
     * @return the minutia rotated around the given center.
     */
    public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {
        //TODO implement

        int x = minutia[1] - centerCol;
        int y = centerRow - minutia[0];

        double rotationRad = Math.toRadians(rotation);

        double newX = x * Math.cos(rotationRad) - y * Math.sin(rotationRad);
        double newY = x * Math.sin(rotationRad) + y * Math.cos(rotationRad);

        int newRow = (int) Math.round(centerRow - newY);
        int newCol = (int) Math.round(newX + centerCol);
        int newOrientation = (minutia[2]+rotation)%360;

        return new int[] {newRow, newCol, newOrientation};
    }

    /**
     * Applies the specified translation to the minutia.
     *
     * @param minutia        the original minutia.
     * @param rowTranslation the translation along the rows.
     * @param colTranslation the translation along the columns.
     * @return the translated minutia.
     */
    public static int[] applyTranslation(int[] minutia, int rowTranslation, int colTranslation) {
        //TODO implement

        minutia[0] = minutia[0] - rowTranslation;
        minutia[1] = minutia[1] - colTranslation;
        // Orientation remains intact

        return minutia;
    }

    /**
     * Computes the row, column, and angle after applying a transformation
     * (translation and rotation).
     *
     * @param minutia        the original minutia.
     * @param centerCol      the column around which the point is rotated.
     * @param centerRow      the row around which the point is rotated.
     * @param rowTranslation the vertical translation.
     * @param colTranslation the horizontal translation.
     * @param rotation       the rotation.
     * @return the transformed minutia.
     */
    public static int[] applyTransformation(int[] minutia, int centerRow, int centerCol, int rowTranslation,
                                            int colTranslation, int rotation) {
        //TODO implement

        // Rotation before Translation
        int[] rotated = applyRotation(minutia, centerRow, centerCol, rotation);

        return applyTranslation(rotated,rowTranslation,colTranslation);
    }

    /**
     * Computes the row, column, and angle after applying a transformation
     * (translation and rotation) for each minutia in the given list.
     *
     * @param minutiae       the list of minutiae.
     * @param centerCol      the column around which the point is rotated.
     * @param centerRow      the row around which the point is rotated.
     * @param rowTranslation the vertical translation.
     * @param colTranslation the horizontal translation.
     * @param rotation       the rotation.
     * @return the list of transformed minutiae.
     */
    public static List<int[]> applyTransformation(List<int[]> minutiae, int centerRow, int centerCol,
                                                  int rowTranslation, int colTranslation, int rotation) {
        //TODO implement

        List<int[]> newMinutiae = new ArrayList<>();

        for (int[] elem : minutiae) {
            int[] newMinutia = applyTransformation(elem, centerRow, centerCol, rowTranslation, colTranslation, rotation);
            newMinutiae.add(newMinutia);
        }

        return newMinutiae;
    }
    /**
     * Counts the number of overlapping minutiae.
     *
     * @param minutiae1      the first set of minutiae.
     * @param minutiae2      the second set of minutiae.
     * @param maxDistance    the maximum distance between two minutiae to consider
     *                       them as overlapping.
     * @param maxOrientation the maximum difference of orientation between two
     *                       minutiae to consider them as overlapping.
     * @return the number of overlapping minutiae.
     */
    public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance,
                                            int maxOrientation) {
        //TODO implement

        // Number of overlaps for the first obtained value
        int overlap = 0;

        for (int[]m1 : minutiae1) {
            for (int[]m2 : minutiae2) {

                int row1 = m1[0];
                int row2 = m2[0];

                int col1 = m1[1];
                int col2 = m2[1];

                int orient1 = m1[2];
                int orient2 = m2[2];

                double sum = Math.pow(row1 - row2, 2) + Math.pow(col1 - col2, 2);
                double Euclid = Math.sqrt(sum);

                int orientDiff = Math.abs(orient1 - orient2);

                if ((Euclid <= maxDistance) && (orientDiff <= maxOrientation)) {
                    overlap++;
                    break; // break allows to avoid unwanted overlaps
                }
            }
        }
        return overlap;
    }

    /**
     * Compares the minutiae from two fingerprints.
     *
     * @param minutiae1 the list of minutiae of the first fingerprint.
     * @param minutiae2 the list of minutiae of the second fingerprint.
     * @return Returns <code>true</code> if they match and <code>false</code>
     *         otherwise.
     */
    public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {
        //TODO implement

        int mFound;
        int maxCnt = 0;

        // Using m2 as rotation center
        for (int[] m1 : minutiae1) {
            for (int[] m2 : minutiae2) {

                int m2Row = m2[0];
                int m2Col = m2[1];
                int rowTranslation = m1[0] - m2[0];
                int colTranslation = m1[1] - m2[1];
                int rotation = m1[2] - m2[2];

                for (int rotOffset = rotation - MATCH_ANGLE_OFFSET
                      ; rotOffset <= rotation + MATCH_ANGLE_OFFSET; rotOffset++) {

                    List<int[]> transformedList = applyTransformation(minutiae1, m2Row, m2Col,
                                                                      rowTranslation, colTranslation, rotOffset);

                    mFound = matchingMinutiaeCount(minutiae2, transformedList,
                                                   DISTANCE_THRESHOLD, ORIENTATION_THRESHOLD);

                    if (mFound >= FOUND_THRESHOLD) {System.out.println("Match Cnt: "+mFound); return true;}
                    else {maxCnt = Math.max(maxCnt,mFound);}

                }
            }
        }
        System.out.println("Match Cnt: "+maxCnt);
        return false;
    }
}