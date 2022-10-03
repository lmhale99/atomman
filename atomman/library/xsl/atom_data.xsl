<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text" omit-xml-declaration="yes" indent="no"/>
  
  <xsl:template match="/*">
    <xsl:apply-templates select="atomic-system"/>
  </xsl:template>

  <xsl:template match="atomic-system">
    
    <!-- # Atoms line -->
    <xsl:value-of select="atoms/natoms"></xsl:value-of>
    <xsl:text> atoms&#xa;</xsl:text>
    
    <!-- # Types line -->
    <xsl:for-each select="atoms/property[name='atype']/data/value">
			<xsl:sort select="." data-type="number" order="descending"/>
			<xsl:if test="position()=1">
				<xsl:value-of select="."/>
			</xsl:if>
		</xsl:for-each>
    <xsl:text> atom types&#xa;</xsl:text>

    <!-- Box dimensions -->
    <xsl:value-of select="format-number(box/origin/value[1], '##0.0000000000000')"/>
    <xsl:text> </xsl:text>
    <xsl:value-of select="format-number(box/origin/value[1] + box/avect/value[1], '##0.0000000000000')"/>
    <xsl:text> xlo xhi&#xa;</xsl:text>

    <xsl:value-of select="format-number(box/origin/value[2], '##0.0000000000000')"/>
    <xsl:text> </xsl:text>
    <xsl:value-of select="format-number(box/origin/value[2] + box/bvect/value[2], '##0.0000000000000')"/>
    <xsl:text> ylo yhi&#xa;</xsl:text>

    <xsl:value-of select="format-number(box/origin/value[3], '##0.0000000000000')"/>
    <xsl:text> </xsl:text>
    <xsl:value-of select="format-number(box/origin/value[3] + box/cvect/value[3], '##0.0000000000000')"/>
    <xsl:text> zlo zhi&#xa;</xsl:text>

    <xsl:if test="box/bvect/value[1]!=0 or box/cvect/value[1]!=0 or box/cvect/value[2]!=0">
      <xsl:value-of select="format-number(box/bvect/value[1], '##0.0000000000000')"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="format-number(box/cvect/value[1], '##0.0000000000000')"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="format-number(box/cvect/value[2], '##0.0000000000000')"/>
      <xsl:text> xy xz yz&#xa;</xsl:text>
    </xsl:if>

    <!-- Atomic values -->
    <xsl:text>&#xa;Atoms # atomic&#xa;&#xa;</xsl:text>
    <xsl:for-each select="atoms/property[name='atype']/data/value">
      <xsl:variable name="row" select="position()"/>
      <xsl:value-of select="$row"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="."/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="format-number(../../../property[name='pos']/data/value[($row - 1) * 3 + 1], '##0.0000000000000')"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="format-number(../../../property[name='pos']/data/value[($row - 1) * 3 + 2], '##0.0000000000000')"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="format-number(../../../property[name='pos']/data/value[($row - 1) * 3 + 3], '##0.0000000000000')"/>
      <xsl:text>&#xa;</xsl:text>
      
    
    </xsl:for-each>

  </xsl:template>
</xsl:stylesheet>